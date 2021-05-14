from os import system
import re
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.type_check import imag
import scipy.linalg, scipy.integrate
from quantum_systems import ODQD, GeneralOrbitalSystem
from IPython.display import display, clear_output

class GHF:
    def __init__(self, l=10, grid_length=10, num_grid_points=201, alpha=1.0, a=0.25, Omega=0.25, omega=2, epsilon0=1, nparticles=2, potential=None, antisymmetrize=True):
        if potential is None:
            self.potential = ODQD.HOPotential(Omega)
        else:
            self.potential = potential
        
        odho = ODQD(l, grid_length, num_grid_points, a, alpha, potential=potential)
        self.system = GeneralOrbitalSystem(n=nparticles, basis_set=odho, anti_symmetrize=antisymmetrize)
        
        self.nparticles = nparticles
        self.omega = omega
        self.epsilon0 = epsilon0
    
    def laser_potential(self, t):
        return self.epsilon0 * np.sin( self.omega * t )
    
    def plot_AO(self): 
        plt.figure(figsize=(8, 5))
        plt.plot(self.system.grid, self.potential(self.system.grid))

        for i in range(int(self.system.l/2)):
            plt.plot( self.system.grid, np.abs(self.system.spf[2*i]) ** 2 + self.system.h[2*i, 2*i].real, label=r"$\chi_{" + f"{i}" + r"}$" )

        plt.grid()
        plt.legend()
        plt.show()
    
    def plot_MO(self, epsilon, C):
        plt.figure(figsize=(8,5))
        plt.plot(self.system.grid, self.potential(self.system.grid))

        for i in range(self.system.l):
            to_plot = np.zeros( self.system.grid.shape, np.complex128 )
            for j in range(self.system.l):
                to_plot += self.system.spf[j] * C[j,i]
            plt.plot( self.system.grid,  np.abs(to_plot)**2 + epsilon[i], label=r"$\phi_{" + f"{i}" + r"}$" )

        plt.grid()
        plt.show()
    
    def plot_overlap(self, time, overlap):
        plt.figure(figsize=(8,5))
        img = plt.imread("theoretical_overlap.png")
        ext = [0.0, 4.0, 0.0, 1.0]
        plt.imshow(img, zorder=0, extent=ext)
        aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
        plt.gca().set_aspect(aspect)
        plt.plot(time, overlap)
    
    def plot_dipole(self, time, dipole):
        plt.figure(figsize=(8,5))
        plt.plot(time, dipole)
        plt.grid()
    
    def plot_energy_per_iteration(self, energy_ghf, energy_rhf):
        plt.figure(figsize=(8,5))
        plt.plot(np.arange(len(energy_ghf)), energy_ghf, label='GHF')
        plt.plot(np.arange(len(energy_rhf)), energy_rhf, label='RHF')
        plt.grid()
        plt.legend()
        
    
    def fill_density_matrix(self, C, full_density=False):
        if full_density is True: 
            imax = len(C)
        else:
            imax = self.nparticles

        density = np.zeros( C.shape, dtype=np.complex128)
        for i in range(imax):
            density += np.outer( np.conjugate(C[:,i]), C[:,i])
        
        return density
    
    def fill_fock_matrix(self, C, t):
        f = np.zeros(self.system.h.shape, dtype=np.complex128)    
        density = self.fill_density_matrix(C)
        
        f = np.einsum('ij,aibj->ab', density, self.system.u, dtype=np.complex128)
        f += self.system.h
        f += self.system.position[0] * self.laser_potential(t)
        return f
    
    def evaluate_total_energy(self, C):
        density = self.fill_density_matrix(C)
        rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
        energy = np.einsum('ij,ij', density, self.system.h) + 0.5*np.einsum('abcd,acbd', rhorho, self.system.u, dtype=np.complex128)
        
        return energy
    
    def eval_one_body_density(self, C, plot_ON=False):
        obd = np.zeros( len(self.system.grid) )
        density = self.fill_density_matrix(C)
        obd = np.einsum('mi,mn,ni->i', self.system.spf, density, self.system.spf, dtype=np.complex128)

        if plot_ON==True: 
            plt.figure(figsize=(10, 10))
            img = plt.imread("theoretical_density.png")
            ext = [-6.0, 6.0, 0.00, 0.4]
            plt.imshow(img, zorder=0, extent=ext)
            aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
            plt.gca().set_aspect(aspect)
            plt.plot(self.system.grid, obd.real)
            plt.grid()
            plt.show()

        return obd
    
    def eval_dipole(self, C):
        density = self.fill_density_matrix(C)
        dipole = np.einsum('mn,mn', density, self.system.position[0])
        return dipole
    
    def solve_TIHF(self, tolerance, max_iter, print_ON=False, energy_per_step_ON=False):
        epsilon, C = scipy.linalg.eigh(np.eye(self.system.l))
        epsilon_old = epsilon
        deltaE=1
        i=0

        if energy_per_step_ON == True:
            energy_per_step = np.zeros((max_iter), dtype=np.complex128)
        else:
            energy_per_step = None
        
        while i<max_iter and deltaE>tolerance:
            f = self.fill_fock_matrix(C, t=0)
            epsilon, C = scipy.linalg.eigh(f)
            deltaE = sum(np.abs(epsilon - epsilon_old))/self.system.l
            epsilon_old = epsilon
            if energy_per_step_ON == True:
                energy_per_step[i] = self.evaluate_total_energy(C)
            i+=1
        
        if deltaE<tolerance and print_ON==True: print("Converged")
        elif i>=max_iter and print_ON==True : print("Max iteration number reached")
        if energy_per_step_ON==True:
            return epsilon, C, energy_per_step
        else:
            return epsilon, C

    def rhsf(self, t, C):
        C = np.reshape(C, self.system.h.shape)
        f = self.fill_fock_matrix(C, t)
        res = -1.0j * (f @ C)
        res = np.reshape(res, (len(res)**2) )
        return res
    
    def solve_TDHF(self, dt, t_max, C0, eval_overlap=False, eval_dipole=False):
        # returns a 1-d array with the right-hand side of time-dependent HF equation
        # this function needs to be defined here because passing f_args to integrator is broken

        i_max = int( np.floor(t_max / dt) )
        time = np.linspace( 0, i_max*dt, i_max ) * 0.5 * self.omega / np.pi

        if eval_overlap==True:
            overlap = np.zeros( (i_max), dtype=np.complex128 )
        if eval_dipole==True:
            dipole = np.zeros( (i_max), dtype=np.complex128 )

        integrator = scipy.integrate.ode(self.rhsf).set_integrator('zvode')
        integrator.set_initial_value( np.reshape(C0, len(C0)**2 ), 0)
                
        i=0
        while integrator.successful() and i<i_max:
            C = integrator.integrate(integrator.t+dt)
            C = np.matrix(np.reshape(C, self.system.h.shape))
            if eval_overlap==True:
                overlap[i] = np.abs( np.linalg.det( C[:,0:2].H @ C0[:,0:2] ))**2
            if eval_dipole==True:
                dipole[i] = self.eval_dipole(C)
            i += 1

        if eval_overlap==True:
            if eval_dipole==True:
                return time, overlap, dipole
            else:
                return time, overlap
        elif eval_overlap==False:
            if eval_dipole==True:
                return time, dipole
            else:
                return time
            
    
    def animation(self, i, line, dt, t_max, C0):
        i_max = int( np.floor(t_max / dt) )
        overlap = np.zeros( (i_max) )
        integrator = scipy.integrate.ode(self.rhsf).set_integrator('zvode')
        integrator.set_initial_value( np.reshape(C0, len(C0)**2 ), 0)

        C = integrator.integrate(integrator.t+dt)
        C = np.matrix(np.reshape(C, self.system.h.shape))
        overlap[i] = np.abs( np.linalg.det( C[:,0:2].H @ C0[:,0:2] ))**2
        
        line.set_ydata(self.eval_one_body_density(C).real)
        return line,
    


class RHF(GHF):
    def __init__(self, l=10, grid_length=10, num_grid_points=201, alpha=1.0, a=0.25, Omega=0.25, omega=2, epsilon0=1, nparticles=2, potential=None, antisymmetrize=False):
        if potential is None:
            self.potential = ODQD.HOPotential(Omega)
        else:
            self.potential = potential
        
        odho = ODQD(l, grid_length, num_grid_points, a, alpha, potential=potential)
        self.system = GeneralOrbitalSystem(n=nparticles, basis_set=odho, anti_symmetrize=antisymmetrize)
        
        self.nparticles = nparticles
        self.omega = omega
        self.epsilon0 = epsilon0
    
    def fill_fock_matrix(self, C, t):
        f = np.zeros(self.system.h.shape, dtype=np.complex128)    
        density = self.fill_density_matrix(C)
        
        f = np.einsum('ij,aibj->ab', density, self.system.u, dtype=np.complex128)
        x = np.zeros(self.system.u.shape)
        x[::2,::2,::2,::2] = 1
        x[1::2,1::2,1::2,1::2] = 1
        x = x*self.system.u
        f -= np.einsum('ij,aijb->ab', density, x, dtype=np.complex128)
        f += self.system.h
        f += self.system.position[0] * self.laser_potential(t)

        return f

    def evaluate_total_energy(self, C):
        density = self.fill_density_matrix(C)
        rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
        energy = np.einsum('ij,ij', density, self.system.h) + 0.5*np.einsum('abcd,acbd', rhorho, self.system.u, dtype=np.complex128)
        x = np.zeros(self.system.u.shape)
        x[::2,::2,::2,::2] = 1
        x[1::2,1::2,1::2,1::2] = 1
        x = x*self.system.u
        energy -= 0.5*np.einsum('abcd,acdb', rhorho, self.system.u, dtype=np.complex128)
        
        return energy
    
    



class wrapfunc(object):
    def __init__(self, func, args=[]):
        self.ff = func
        self.args = args
    def f(self, t, C):
        return self.ff(t, C, *self.args)

    