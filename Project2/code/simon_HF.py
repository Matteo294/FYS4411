import numpy as np
import numba
import matplotlib.pyplot as plt
from quantum_systems import ODQD, GeneralOrbitalSystem
import scipy
def construct_Density_matrix(C,number_electrons,l):
    n=np.eye(2*l)
    for i in range(number_electrons,2*l):
        n[i,i]=0
    P=C@n@C.conjugate().T
    #return P
    """
    P=np.zeros((2*l,2*l),dtype=np.complex64)
    for tau in range(2*l):
        for sigma in range(2*l):
            for i in range(number_electrons):
                P[tau,sigma]=np.conjugate(C[sigma,i])*C[tau,i]
    """
    return P

def costruct_Fock_matrix(P,l,number_electrons,system,anti_symmetrize=True):
    udirect=np.zeros((2*l,2*l),dtype=np.complex64)
    uexchange=np.zeros((2*l,2*l),dtype=np.complex64)
    for tau in range(2*l):
        for sigma in range(2*l):
            for mu in range(2*l):
                for nu in range(2*l):
                    udirect[mu,nu]+=P[tau,sigma]*system.u[mu,sigma,nu,tau]
                    if not anti_symmetrize:
                        uexchange[mu,nu]+=P[tau,sigma]*system.u[mu,sigma,tau,nu]
    F=system.h+udirect-uexchange
    return F
def solve(system,number_electrons,number_basissets,C=None,anti_symmetrize=True,tol=1e-8,maxiter=100):
    if C is None:
        C=np.random.rand(2*number_basissets,2*number_basissets)
    P=np.zeros(C.shape)
    converged=False
    for i in range(maxiter):
        if(i>=1):
            convergence_difference=np.max(np.abs(P-P_old))
            if (convergence_difference<tol):
                converged=True
                break
        P_old=P
        P=construct_Density_matrix(C,number_electrons,number_basissets)
        print(P[0,0])
        F=costruct_Fock_matrix(P,number_basissets,number_electrons,system,anti_symmetrize)
        epsilon, C = scipy.linalg.eigh(F)
        #print(epsilon[0])
        #print(C[0,0:2])
    return F,epsilon, C, converged
def compute_energy(C,F,system,number_electrons,number_basissets):
    energy=0
    P=construct_Density_matrix(C,number_electrons,number_basissets)
    for mu in range (2*number_basissets):
        for nu in range(2*number_basissets):
            energy+=P[nu,mu]*(system.h[nu,mu]+F[nu,mu])
    return energy*0.5


l=10
grids_length=10
num_grid_points=101
omega=0.25
a=0.25
odho = ODQD(l, grids_length, num_grid_points, a=a, alpha=1, potential=ODQD.HOPotential(omega))
number_electrons=2
anti_symmetrize=True
system=GeneralOrbitalSystem(n=number_electrons,basis_set=odho,anti_symmetrize=anti_symmetrize)
print("Reference energy: ",(system.compute_reference_energy()))
#epsilon,C=scipy.linalg.eigh(system.construct_fock_matrix(system.h,system.u)); #Øyvinds kode
#print(epsilon) #epsilon fra Øyvinds kode
C=np.eye(2*l) #Create an initial guess for the coefficients
F,epsilon,C,converged=solve(system,number_electrons,l,anti_symmetrize=anti_symmetrize,tol=1e-4,maxiter=10,C=C)
#print(np.diag(F)-epsilon) # This should be as close to zero as possible
#print("Converged: ",converged)
#print("Energy: ", compute_energy(C,F,system,number_electrons,l))
print(epsilon)
