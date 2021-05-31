from logging import log
from operator import pos
from os import system
import re
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import shape
from numpy.lib.type_check import imag
import scipy.linalg, scipy.integrate
from quantum_systems import ODQD, GeneralOrbitalSystem
from IPython.display import display, clear_output
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

class GHF:
    
    def __init__(self, l=10, grid_length=10, num_grid_points=201, alpha=1.0, a=0.25, Omega=0.25, omega=2, epsilon0=1, nparticles=2, potential=None, antisymmetrize=True):
        if potential is None:
            self.potential = ODQD.HOPotential(Omega)
        else:
            self.potential = potential
        
        odho = ODQD(l, grid_length, num_grid_points, a, alpha, potential=self.potential)
        self.system = GeneralOrbitalSystem(n=nparticles, basis_set=odho, anti_symmetrize=antisymmetrize)
        
        self.nparticles = nparticles
        self.Omega = Omega
        self.omega = omega
        self.epsilon0 = epsilon0
        
    def laser_potential(self, t):
        return self.epsilon0 * np.sin( self.omega * t )

    def fill_density_matrix(self, C):
        density = np.zeros( C.shape, dtype=np.complex128)
        for i in range(self.nparticles):
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
            plt.figure(figsize=(13, 13))
            img = plt.imread("theoretical_density.png")
            ext = [-6.0, 6.0, 0.00, 0.4]
            plt.imshow(img, zorder=0, extent=ext)
            aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
            plt.gca().set_aspect(aspect)
            plt.plot(self.system.grid, obd.real)
            plt.xlabel('x', fontsize=16)
            plt.ylabel(r'$\rho(x)$', fontsize=16)
            plt.xlim([-10, 10])
            ax = plt.gca()
            ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
            plt.grid()
            plt.show()

        return obd
    
    def eval_dipole(self, C):
        density = self.fill_density_matrix(C)
        dipole = np.einsum('mn,mn', density, self.system.position[0])
        return dipole
    
    def solve_TIHF(self, tolerance, max_iter, print_ON=False, eval_energy_per_step=False, eval_delta=False):
        epsilon, C = scipy.linalg.eigh(np.eye(self.system.l))
        epsilon_old = epsilon
        deltaE=1
        i=0

        if eval_energy_per_step == True:
            energy_per_step = np.zeros((max_iter), dtype=np.complex128)
        else:
            energy_per_step=None
        
        if eval_delta==True:
            delta_per_step = np.zeros((max_iter), dtype=np.complex128)
        else:
            delta_per_step=None
        
        while i<max_iter and deltaE>tolerance:
            f = self.fill_fock_matrix(C, t=0)
            epsilon, C = scipy.linalg.eigh(f)
            deltaE = np.sum(np.abs(epsilon - epsilon_old))/self.system.l
            epsilon_old = epsilon
            if eval_energy_per_step == True:
                energy_per_step[i] = self.evaluate_total_energy(C)
            if eval_delta==True:
                delta_per_step[i] = deltaE
            i+=1
        
        if deltaE<tolerance and print_ON==True: print("Converged")
        elif i>=max_iter and print_ON==True : print("Max iteration number reached")
        
        return epsilon, C, energy_per_step, delta_per_step
        
            

    def rhsf(self, t, C):
        C = np.reshape(C, self.system.h.shape)
        f = self.fill_fock_matrix(C, t)
        res = -1.0j * (f @ C)
        res = np.reshape(res, (len(res)**2) )
        return res
    
    def rhsf_OFF(self, t, C):
        C = np.reshape(C, self.system.h.shape)
        f = self.fill_fock_matrix(C, t=0)
        res = -1.0j * (f @ C)
        res = np.reshape(res, (len(res)**2) )
        return res
    
    def solve_TDHF(self, tstart, dt, t_max, C0, eval_overlap=False, eval_dipole=False, laser_ON=True):
        i_max = int( np.floor((t_max - tstart)/ dt) )
        time = np.linspace( tstart, tstart+i_max*dt, i_max )

        if eval_overlap==True:
            overlap = np.zeros( (i_max), dtype=np.complex128 )
        else: 
            overlap=None
        
        if eval_dipole==True:
            dipole = np.zeros( (i_max), dtype=np.complex128 )
        else:
            dipole=None

        if laser_ON==True:
            integrator = scipy.integrate.ode(self.rhsf).set_integrator('zvode')
        elif laser_ON==False:
            integrator = scipy.integrate.ode(self.rhsf_OFF).set_integrator('zvode')
        
        integrator.set_initial_value( np.reshape(C0, (len(C0)**2,) ), tstart)

        i=0
        while integrator.successful() and i<i_max:
            C = integrator.integrate(integrator.t+dt)
            C = np.matrix(np.reshape(C, self.system.h.shape))

            if eval_overlap==True:
                overlap[i] = np.abs(np.linalg.det( C[:,0:2].H @ C0[:,0:2] ))**2
            if eval_dipole==True:
                dipole[i] = self.eval_dipole(C)
            i += 1
            if i%10==0:
                print('\rprogress = %.2f' % (100*i/i_max), end='')
        
        print('\r', end='')
        C = np.array(C)

        return C, time, overlap, dipole

    def fourier_analysis(self, tolerance, max_iter, t_laser_ON, t_max, dt):
        epsilon, C0, energy_per_step, delta_per_step = self.solve_TIHF(tolerance=tolerance, max_iter=max_iter, print_ON=False)
        C1, time1, overlap1, dipole1 = self.solve_TDHF(0, dt, t_laser_ON, C0, eval_overlap=True, eval_dipole=True, laser_ON=True)
        C2, time2, overlap2, dipole2 = self.solve_TDHF(t_laser_ON,  dt, t_max, C1, eval_overlap=True, eval_dipole=True, laser_ON=False)
        time = np.concatenate((time1,time2))
        dipole = np.concatenate((dipole1, dipole2))
        overlap = np.concatenate((overlap1, overlap2))

        dipoleFFT = np.fft.fft(dipole2)
        dipolefreqFFT = np.fft.fftfreq(len(dipole2), dt)
        overlapFFT = np.fft.fft(overlap2)
        overlapfreqFFT = np.fft.fftfreq(len(overlap2), dt)       
        
        return C2, time, dipole, overlap, dipoleFFT, dipolefreqFFT, overlapFFT, overlapfreqFFT


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
        
        self.system = ODQD(l, grid_length, num_grid_points, a, alpha, potential=self.potential)
        self.nparticles = nparticles
        self.Omega = Omega
        self.omega = omega
        self.epsilon0 = epsilon0
    
    def fill_density_matrix(self, C, full_density=False):
        if full_density is True: 
            imax = len(C)
        else:
            imax = int(self.nparticles/2)

        density = np.zeros( C.shape, dtype=np.complex128)
        for i in range(imax):
            density += np.outer( np.conjugate(C[:,i]), C[:,i])
        
        return density
        return 

    def fill_fock_matrix(self, C, t):
        f = np.zeros(self.system.h.shape, dtype=np.complex128)    
        density = self.fill_density_matrix(C)
        f = 2*np.einsum('ij,aibj->ab', density, self.system.u, dtype=np.complex128)
        f -= np.einsum('ij,aijb->ab', density, self.system.u, dtype=np.complex128)
        f += self.system.h
        f += self.system.position[0] * self.laser_potential(t)

        return f
    
    def evaluate_total_energy(self, C):
        density = self.fill_density_matrix(C)
        rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
        energy = 2*np.einsum('ij,ij', density, self.system.h) 
        energy += 2*np.einsum('abcd,acbd', rhorho, self.system.u, dtype=np.complex128)
        energy -= np.einsum('abcd,acdb', rhorho, self.system.u, dtype=np.complex128)
        
        return energy
    
    def eval_one_body_density(self, C, plot_ON=False):
        obd = np.zeros( len(self.system.grid) )
        density = self.fill_density_matrix(C)
        obd = 2*np.einsum('mi,mn,ni->i', self.system.spf, density, self.system.spf, dtype=np.complex128)

        if plot_ON==True: 
            plt.figure(figsize=(13, 13))
            img = plt.imread("theoretical_density.png")
            ext = [-6.0, 6.0, 0.00, 0.4]
            plt.imshow(img, zorder=0, extent=ext)
            aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
            plt.gca().set_aspect(aspect)
            plt.plot(self.system.grid, obd.real)
            plt.grid()
            plt.show()

        return obd

    