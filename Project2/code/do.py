import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg, scipy.integrate
from quantum_systems import ODQD, GeneralOrbitalSystem
from IPython.display import display, clear_output

#can be called only after creation of odho and before creation of system!!!
# plots the Atormic Orbital basis --> in our case the HO eigenstates
def plot_AO(odho, potential): 
    plt.figure(figsize=(8, 5))
    plt.plot(odho.grid, potential(odho.grid))

    for i in range(odho.l):
        plt.plot( odho.grid, np.abs(odho.spf[i]) ** 2 + odho.h[i, i].real, label=r"$\chi_{" + f"{i}" + r"}$" )

    plt.grid()
    plt.legend()
    plt.show()

# prints the Molecular Orbitals obtained after the convergence of the HF method
# the basis has not been changed for the system
# it is used just for a comparison with the output of plot_MO_abc()
def plot_MO(system, potential, epsilon, C):
    plt.figure(figsize=(8,5))
    plt.plot(system.grid, potential(system.grid))

    for i in range(system.l):
        to_plot = np.zeros( system.grid.shape, np.complex128 )
        for j in range(system.l):
            to_plot += system.spf[j] * C[j,i]
        plt.plot( system.grid,  np.abs(to_plot)**2 + epsilon[i], label=r"$\phi_{" + f"{i}" + r"}$" )

    plt.grid()
    plt.legend()
    plt.show()

# abc = after basis change, can be called only after system.change_basis(C)
# plots the Molecular Orbitals obtained after the convergence of the HF procedure
# here the basis of the system has already been changed
def plot_MO_abc(system, potential, epsilon):
    plt.figure(figsize=(8, 5))
    plt.plot(system.grid, potential(system.grid))

    for i in range(system.l):
        plt.plot( system.grid, np.abs(system.spf[i]) ** 2 + epsilon[i], label=r"$\phi_{" + f"{i}" + r"}$" )

    plt.grid()
    plt.legend()
    plt.show()

# changes basis and modifies system as a consequence
def change_basis(system, C):
    
    ####### system.spf ########
    spf = np.copy(system.spf)
    for i in range(system.l):
        system.spf[i] = np.zeros(system.grid.shape, dtype=np.complex128)
        for j in range(system.l):
            system.spf[i] += spf[j,:] * C[j,i]   

def laser_potential(t, omega, epsilon0):
    return epsilon0 * np.sin( omega * t )

def fill_fock_matrix(system, nparticles, C, t, omega, epsilon0):
    f = np.zeros(system.h.shape, dtype=np.complex128)    
    density = fill_density_matrix(nparticles, C)
    
    f = np.einsum('ij,aibj->ab', density, system.u, dtype=np.complex128)
    f += system.h
    f += system.position[0] * laser_potential(t,omega,epsilon0)
    return f


# density matrix
def fill_density_matrix(nparticles, coefficients):
    density = np.zeros( coefficients.shape, dtype=np.complex128)
    for i in range(nparticles):
        density += np.outer( np.conjugate(coefficients[:,i]), coefficients[:,i])
    
    return density


# contains two THEORETICALLY equivalent ways of evaluating the energy of the system
# at the moment it works only after the basis has been changed
def evaluate_total_energy(system, epsilon=None):
    energy = 0
    
    if epsilon is None:
        for i in range(2):
            energy += system.h[i,i]
            for j in range(2):
                energy += 0.5 * system.u[i,j,i,j]
    else:
        for i in range(2):
            for j in range(2):
                energy -= system.u[i,j,i,j] 
        energy = energy*0.5 + epsilon[0] + epsilon[1]

    return energy


def eval_one_body_density(system, nparticles, coefficients):
    obd = np.zeros( len(system.grid) )
    density = fill_density_matrix(nparticles, coefficients)
    obd = np.einsum('mi,mn,ni->i', system.spf, density, system.spf, dtype=np.complex128)
    return obd


# evaluates the one-body density and plots it overlapped with the curve from the article
def plot_overlap_one_body_density(system, obd):  
    plt.figure(figsize=(16, 16))
    img = plt.imread("theoretical_density.png")
    ext = [-6.0, 6.0, 0.00, 0.4]
    plt.imshow(img, zorder=0, extent=ext)
    aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
    plt.gca().set_aspect(aspect)
    plt.plot(system.grid, obd.real)
    plt.grid()
    plt.show()

def plot_overlap_slater_det(system, overlap, time):
    plt.figure(figsize=(12,8))
    img = plt.imread("theoretical_overlap.png")
    ext = [0.0, 4.0, 0.0, 1.0]
    plt.imshow(img, zorder=0, extent=ext)
    aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
    plt.gca().set_aspect(aspect)
    plt.plot(time, overlap)

def solve_TIHF(system, nparticles, tolerance, max_iter, print_on=False):
    epsilon, C = scipy.linalg.eigh(system.h)
    epsilon_old = epsilon
    deltaE=1
    i = 0
    
    while i<max_iter and deltaE>tolerance:
        f = fill_fock_matrix(system, nparticles, C, t=0, omega=0, epsilon0=0)
        epsilon, C = scipy.linalg.eigh(f)
        deltaE = sum(np.abs(epsilon - epsilon_old))/system.l
        epsilon_old = epsilon
        i+=1
    
    if deltaE<tolerance and print_on==True: print("Converged")
    elif i>=max_iter and print_on==True : print("Max iteration number reached")
    return epsilon, C


def solve_TDHF(system, dt, t_max, C0, omega, epsilon0, nparticles, animation=False):
    # returns a 1-d array with the right-hand side of time-dependent HF equation
    # this function needs to be defined here because passing f_args to integrator is broken
    def rhs(t, C):
        C = np.reshape(C, system.h.shape)
        f = fill_fock_matrix(system, nparticles, C, t, omega, epsilon0)
        res = -1.0j * (f @ C)
        res = np.reshape(res, (len(res)**2) )
        return res

    i_max = int( np.floor(t_max / dt) )
    overlap = np.zeros( (i_max) )
    time = np.linspace( 0, i_max*dt, i_max ) * 0.5 * omega / np.pi

    integrator = scipy.integrate.ode(rhs).set_integrator('zvode')
    integrator.set_initial_value( np.reshape(C0, len(C0)**2 ), 0)
    obd0 = eval_one_body_density(system, nparticles, C0).real

    if animation==True:
        fig = plt.figure(figsize=(8,5))
        ax = fig.add_subplot(1, 1, 1)   
    
    
    i=0
    while integrator.successful() and i<i_max:
        C = integrator.integrate(integrator.t+dt)
        C = np.matrix(np.reshape(C, system.h.shape))
        overlap[i] = np.abs( np.linalg.det( C[:,0:2].H @ C0[:,0:2] ))**2
        if animation==True and i%(int(i_max*0.01))==0:
            ax.cla()
            ax.plot(system.grid, obd0, color='red')
            ax.plot(system.grid, eval_one_body_density(system, nparticles, C).real)
            display(fig)
            clear_output(wait=True)
        i += 1
    
    return overlap, time
    