import MPSPyLib as mps
import numpy as np
from sys import version_info
import matplotlib.pyplot as plt
import sys

 # Build operators
Operators = mps.BuildBoseOperators(5)
Operators['interaction'] = 0.5 * (np.dot(Operators['nbtotal'], Operators['nbtotal']) - Operators['nbtotal'])

# Define Hamiltonian MPO
H = mps.MPO(Operators)
H.AddMPOTerm('bond', ['bdagger', 'b'], hparam='t', weight=-1.0)
H.AddMPOTerm('site', 'interaction', hparam='U', weight=1.0)

myObservables = mps.Observables(Operators)
myObservables.AddObservable('site', 'nbtotal', 'n')
myObservables.AddObservable('Lambda', [])

myConv = mps.MPSConvParam(max_bond_dimension=20, max_num_sweeps=6)
myConv.AddModifiedConvergenceParameters(0, ['max_bond_dimension',
                                                'local_tol'], [50, 1E-14])

myKrylovConv = mps.KrylovConvParam(max_num_lanczos_iter=20,lanczos_tol=1E-6)

dynObservables = mps.Observables(Operators)
dynObservables.AddObservable('site', 'nbtotal', name='n')
dynObservables.AddObservable('Lambda', [])

tau=10
U = 10.0
t = 1.0
Parameters = []
L = 6


Quenches = mps.QuenchList(H)

def Ufuncdown(t):
     return U

Quenches.AddQuench(['U'], tau, min(tau / 100.0, 0.1),
                           [Ufuncdown], ConvergenceParameters=myKrylovConv)

Parameters.append({
            'simtype'                   : 'Finite',
            # Directories
            'job_ID'                    : 'Bose_Hubbard_mod',
            'unique_ID'                 : 'tau_'+ str(tau),
            'Write_Directory'           : 'TMP_03/',
            'Output_Directory'          : 'OUTPUTS_03/',
            # System size and Hamiltonian parameters
            'L'                         : L,
            't'                         : t, 
            'U'                         : U, 
            # Specification of symmetries and good quantum numbers
            'Abelian_generators'        : ['nbtotal'],
            'Abelian_quantum_numbers'   : [L],
            # Convergence parameters
            'verbose'                   : 1,
            'logfile'                   : True,
            'MPSObservables'            : myObservables,
            'MPSConvergenceParameters'  : myConv,
            'Quenches'                  : Quenches,
            'DynamicsObservables'       : dynObservables
})

PostProcess='T'
MainFiles = mps.WriteFiles(Parameters, Operators, H,PostProcess=PostProcess)
if(not PostProcess):
        if os.path.isfile('./Execute_MPSMain'):
            RunDir = './'
        else:
            RunDir = None
        mps.runMPS(MainFiles, RunDir=RunDir)
		#return
Outputs = mps.ReadStaticObservables(Parameters)
DynOutputs = mps.ReadDynamicObservables(Parameters)

print(DynOutputs[0][0]['bond_dimension'])
