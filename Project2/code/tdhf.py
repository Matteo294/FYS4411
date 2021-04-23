import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.integrate
from quantum_systems import ODQD, GeneralOrbitalSystem

def plot_AO(odho):
    fig = plt.figure(figsize=(8, 5))
    plt.plot(odho.grid, potential(odho.grid))

    for i in range(l):
        plt.plot( odho.grid, np.abs(odho.spf[i]) ** 2 + odho.h[i, i].real, label=r"$\chi_{" + f"{i}" + r"}$" )

    plt.grid()
    plt.legend()
    plt.show()

def solve_HF(system, tolerance, max_iter):
    i = 0
    f = np.zeros( system.h.shape, dtype=np.complex128 )
    epsilon, C = scipy.linalg.eigh(system.h)
    print(C[0,0])
    deltaE=1
    
    while i<max_iter and deltaE>tolerance:
        
        for j in range(len(system.h)):
            for k in range(len(system.h)):
                
                two_body_term = 0
                for m in range(len(system.h)):
                    for n in range(len(system.h)): 
                        if (m+n)%2!=0:
                            two_body_term += 0
                        else:
                            sum_c = 0
                            for p in range(len(system.h)):
                                sum_c += C[m,p]*C[n,p]
                            two_body_term +=  sum_c*system.u[j,m,k,n]
                f[j,k] = two_body_term + system.h[j,k]
                
        epsilon, C = scipy.linalg.eigh(f)
        print(epsilon)
        i+=1


                            
### PARAMETERS ###
l = 5 # Number of eigenstates of the HO potential --> we use these functions to generate the single particle WF
grid_length = 10  # The width of the one-dimensional grid
num_grid_points = 101  # The number of discretized points on the grid.
# More points give better results for the single-particle basis at the cost of slower setup.
alpha = 1  # The strength of the Coulomb interaction 
a = 0.25  # The shielding parameter in the Coulomb interaction potential
omega = 0.25  # The frequency of the harmonic oscillator trap
potential=ODQD.HOPotential(omega)
nparticles = 2

odho = ODQD(l, grid_length, num_grid_points, a=a, alpha=alpha, potential=potential )
system = GeneralOrbitalSystem(n=nparticles, basis_set=odho, anti_symmetrize=True)
solve_HF(system, tolerance=1e-3, max_iter=10)




