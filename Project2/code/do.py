import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg, scipy.integrate
from quantum_systems import ODQD, GeneralOrbitalSystem
from IPython.display import display, clear_output


def plot_AO(odho, potential): 

    """ This function can be called only after creation of odho and before creation of system.
    This plots the Atormic Orbital basis --> in our case the HO eigenstates """


    plt.figure(figsize=(8, 5))
    plt.plot(odho.grid, potential(odho.grid))

    for i in range(odho.l):
        plt.plot( odho.grid, np.abs(odho.spf[i]) ** 2 + odho.h[i, i].real, label=r"$\chi_{" + f"{i}" + r"}$" )

    plt.grid()
    plt.legend()
    plt.show()


def plot_MO(system, potential, epsilon, C):


    """ This function prints the Molecular Orbitals obtained after the convergence of the HF method
    the basis has not been changed for the system it is used just for a comparison with the output of plot_MO_abc()"""


    plt.figure(figsize=(8,5))
    plt.plot(system.grid, potential(system.grid))

    for i in range(system.l):
        to_plot = np.zeros( system.grid.shape, np.complex128 )
        for j in range(system.l):
            to_plot += system.spf[j] * C[j,i]
        plt.plot( system.grid,  np.abs(to_plot)**2 + epsilon[i], label=r"$\phi_{" + f"{i}" + r"}$" )

    plt.grid()
    #plt.legend()
    plt.show()

def plot_MO_abc(system, potential, epsilon):
    """abc = after basis change, can be called only after system.change_basis(C).
    This function plots the Molecular Orbitals obtained after the convergence of the HF procedure when the basis of the system has already been changed"""
    plt.figure(figsize=(8, 5))
    plt.plot(system.grid, potential(system.grid))

    for i in range(system.l):
        plt.plot( system.grid, np.abs(system.spf[i]) ** 2 + epsilon[i], label=r"$\phi_{" + f"{i}" + r"}$" )

    plt.grid()
    #plt.legend()
    plt.show()

def change_basis(system, C):


    """changes basis and modifies system as a consequence (change basis of system.spf , system.h, system.u)
    """


    ####### system.spf ########
    #spf = np.copy(system.spf)
    spf = system.spf
    system._basis_set._spf = np.einsum('ki,kj->ij', C, spf)
    
    ####### system.h #########
    h = np.copy(system.h)
    cc = np.einsum('ai,bj->aibj', np.conjugate(C), C)
    system._basis_set._h = np.einsum('aibj,ab->ij', cc, h)

    ####### system.u #########
    u = np.copy(system.u)
    u1 = np.einsum('am,abcd->mbcd', np.conjugate(C), u)
    u2 = np.einsum('bn, mbcd->mncd', np.conjugate(C), u1)
    u3 = np.einsum('cp, mncd->mnpd', C, u2)
    system._basis_set._u = np.einsum('dq,mnpd->mnpq', C, u3)

def laser_potential(t, omega, epsilon0):
    """ This function sets up the laser potential"""
    return epsilon0 * np.sin( omega * t )


# density matrix
def fill_density_matrix(coefficients, Nparticles=None):
    """This function fills the density matrix from the coefficient matrix""" 
    if Nparticles is None: 
        Nparticles=len(coefficients)

    density = np.zeros( coefficients.shape, dtype=np.complex128)
    for i in range(Nparticles):
        density += np.outer( np.conjugate(coefficients[:,i]), coefficients[:,i])
    
    return density


def fill_fock_matrix(system, nparticles, C, t, omega, epsilon0, unrestricted=False):
    """This function fills the fock matrix""" 
    f = np.zeros(system.h.shape, dtype=np.complex128)    
    density = fill_density_matrix(C,nparticles)
    
    if unrestricted==False:
        f = np.einsum('ij,aibj->ab', density, system.u, dtype=np.complex128)
        f += system.h
        f += system.position[0] * laser_potential(t, omega, epsilon0)
    else:
        f = np.einsum('ij,aibj->ab', density, system.u, dtype=np.complex128)
        x = np.zeros(system.u.shape)
        x[::2,::2,::2,::2] = 1
        x[1::2,1::2,1::2,1::2] = 1
        x = x*system.u
        f -= np.einsum('ij,aijb->ab', density, x, dtype=np.complex128)
        f += system.h
        f += system.position[0] * laser_potential(t, omega, epsilon0)

    return f


# evaluates the total energy of the system
# adapts the formulas depending on the fact that the basis change has been performed or not
def evaluate_total_energy(system, nparticles, C, changed=False, epsilon=None, unrestricted=False):
    """This function evaluates the total energy of the system. It adapts the formulas depending on the fact that the basis change has been performed or not"""
    energy = 0
    if unrestricted==False:
        if changed is True:
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
        # HO basis
        else:
            density = fill_density_matrix(C,nparticles)
            rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
            energy = np.einsum('ij,ij', density, system.h) + 0.5*np.einsum('abcd,acbd', rhorho, system.u, dtype=np.complex128)
    
    else:
        density = fill_density_matrix(C,nparticles)
        rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
        energy = np.einsum('ij,ij', density, system.h) + 0.5*np.einsum('abcd,acbd', rhorho, system.u, dtype=np.complex128)
        x = np.zeros(system.u.shape)
        x[::2,::2,::2,::2] = 1
        x[1::2,1::2,1::2,1::2] = 1
        x = x*system.u
        energy -= 0.5*np.einsum('abcd,acdb', rhorho, system.u, dtype=np.complex128)

    
    return energy


def eval_one_body_density(system, nparticles, coefficients):
    """This function evvaluates the body density """
    obd = np.zeros( len(system.grid) )
    density = fill_density_matrix(coefficients,nparticles)
    obd = np.einsum('mi,mn,ni->i', system.spf, density, system.spf, dtype=np.complex128)
    return obd


# 
def plot_overlap_one_body_density(system, obd):  
    """This function plots the one body density overlapped with the curve from the article (put link)"""
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

    """ This function plots the one body density overlapped with the curve from the article (put link)"""

    plt.figure(figsize=(12,8))
    img = plt.imread("theoretical_overlap.png")
    ext = [0.0, 4.0, 0.0, 1.0]
    plt.imshow(img, zorder=0, extent=ext)
    aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
    plt.gca().set_aspect(aspect)
    plt.plot(time, overlap)

def solve_TIHF(system, nparticles, tolerance, max_iter, print_on=False, unrestricted=False):
    """ This function solves the time Independent hartree fock"""
    epsilon, C = scipy.linalg.eigh(system.h)
    epsilon_old = epsilon
    deltaE=1
    i = 0
    energy = np.zeros((max_iter), dtype=np.complex128)
    
    while i<max_iter and deltaE>tolerance:
        f = fill_fock_matrix(system, nparticles, C, t=0, omega=0, epsilon0=0, unrestricted=unrestricted)
        epsilon, C = scipy.linalg.eigh(f)
        deltaE = sum(np.abs(epsilon - epsilon_old))/system.l
        epsilon_old = epsilon
        energy[i] = evaluate_total_energy(system, nparticles, C, changed=False, epsilon=None, unrestricted=unrestricted)
        i+=1
    
    if deltaE<tolerance and print_on==True: print("Converged")
    elif i>=max_iter and print_on==True : print("Max iteration number reached")
    return epsilon, C, energy

class wrapfunc(object):
    """ wrapfunction"""
    def __init__(self, func, args=[]):
        self.ff = func
        self.args = args
    def f(self, t, C):
        return self.ff(t, C, *self.args)

def rhsf(t, C, system, omega, epsilon0, nparticles):
    """Right half side"""
    C = np.reshape(C, system.h.shape)
    f = fill_fock_matrix(system, nparticles, C, t, omega, epsilon0)
    res = -1.0j * (f @ C)
    res = np.reshape(res, (len(res)**2) )
    return res

def solve_TDHF(system, dt, t_max, C0, omega, epsilon0, nparticles, integrator, animation=False):
    """ This function solves the time dependent hartree fock and 
    returns a 1-d array with the right-hand side of time-dependent HF equation.
    This function needs to be defined here because passing f_args to integrator is broken"""

    i_max = int( np.floor(t_max / dt) )
    overlap = np.zeros( (i_max) )
    time = np.linspace( 0, i_max*dt, i_max ) * 0.5 * omega / np.pi
    
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

def animation(i, line, system, dt, t_max, C0, omega, epsilon0, nparticles, integrator):
    """This function plots the animation of the one body density over time"""
    i_max = int( np.floor(t_max / dt) )
    overlap = np.zeros( (i_max) )
    time = np.linspace( 0, i_max*dt, i_max ) * 0.5 * omega / np.pi

    obd0 = eval_one_body_density(system, nparticles, C0).real

    C = integrator.integrate(integrator.t+dt)
    C = np.matrix(np.reshape(C, system.h.shape))
    overlap[i] = np.abs( np.linalg.det( C[:,0:2].H @ C0[:,0:2] ))**2
    
    line.set_ydata(eval_one_body_density(system, nparticles, C).real)
    return line,
    