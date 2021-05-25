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

    def plot_onebody_comparison(self, onebodyg, onebodyr, save_fig=False):
        plt.figure(figsize=(16, 16))
        img = plt.imread("theoretical_density.png")
        ext = [-6.0, 6.0, 0.00, 0.4]
        plt.imshow(img, zorder=0, extent=ext)
        aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
        plt.gca().set_aspect(aspect)
        plt.plot(self.system.grid, onebodyg.real, label='GHF')
        plt.plot(self.system.grid, onebodyr.real, label='RHF')
        plt.xlabel('x', fontsize=16)
        plt.ylabel(r'$\rho(x)$', fontsize=16)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
        plt.grid()
        plt.legend(fontsize=12)
        if save_fig==True:
            plt.savefig('../paper/images/onebody_density_comp_article.pdf', bbox_inches='tight')
        plt.show()
    
    def plot_energy_per_iteration(self, energy_ghf, energy_rhf, save_fig=False):
        plt.figure(figsize=(10,7))
        plt.plot(np.arange(len(energy_ghf)), energy_ghf.real, label='GHF')
        plt.plot(np.arange(len(energy_rhf)), energy_rhf.real, label='RHF')
        plt.xlabel('Iteration', fontsize=22)
        plt.ylabel('Energy', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        plt.legend(fontsize=14)

        if save_fig==True:
            plt.savefig('../paper/images/energy_at_every_step.pdf', bbox_inches='tight')
        plt.show()
        
    def plot_delta_per_iteration(self, delta_ghf, delta_rhf, tolerance, save_fig=False):
        plt.figure(figsize=(10,7))
        plt.plot(np.arange(len(delta_ghf)), delta_ghf.real, label='GHF')
        plt.plot(np.arange(len(delta_rhf)), delta_rhf.real, label='RHF')
        plt.plot(np.arange(len(delta_rhf)), tolerance*np.ones(len(delta_ghf)), linestyle='--', color='red', label=r'$\delta$')
        plt.xlabel('Iteration', fontsize=22)
        plt.ylabel(r'$\Delta$', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        ax.set_yscale('log')
        plt.grid()
        plt.legend(fontsize=14)

        if save_fig==True:
            plt.savefig('../paper/images/delta_at_every_step.pdf', bbox_inches='tight')
        plt.show()
    
    def plot_overlap(self, time, overlap, save_fig=False):
        plt.figure(figsize=(16,5))
        img = plt.imread("theoretical_overlap.png")
        ext = [0.0, 4.0, 0.0, 1.0]
        plt.imshow(img, zorder=0, extent=ext)
        aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
        plt.gca().set_aspect(aspect)
        plt.plot(time* 0.5 * self.omega / np.pi, overlap)
        plt.xlabel(r'$t\omega/2\pi}$', fontsize=16)
        plt.ylabel(r'$\xi(t)$', fontsize=16)
        #plt.xlim([0,8])
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
        plt.grid()
        if save_fig==True:
            plt.savefig('../paper/images/overlap_comp_article.pdf', bbox_inches='tight')
        plt.show()

    def plot_dipole(self, time, dipole, save_fig=False):
        plt.figure(figsize=(10,7))
        plt.plot(time* 0.5 * self.omega / np.pi, dipole)
        plt.xlabel(r'$t\omega/2\pi}$', fontsize=22)
        plt.ylabel(r'$\overline{x}(t)$', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        if save_fig==True:
            plt.savefig('../paper/images/time_dependent_dipole.pdf', bbox_inches='tight')
        plt.show()
        
    def plot_overlaps(self, time, overlap1, overlap2, overlap3, save_fig=False):
        plt.figure(figsize=(16,5))
        plt.plot(time* 0.5 * self.Omega / np.pi, overlap1, label=r'$\omega=8\Omega$')
        plt.plot(time* 0.5 * self.Omega / np.pi, overlap2, label=r'$\omega=4\Omega$')
        plt.plot(time* 0.5 * self.Omega / np.pi, overlap3, label=r'$\omega=32\Omega$')
        plt.xlabel(r'$t\Omega/2\pi$', fontsize=16)
        plt.ylabel(r'$\xi(t)$', fontsize=16)
        plt.xlim([0, 1.2])
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
        plt.grid()
        plt.legend(fontsize=16)
        if save_fig==True:
            plt.savefig('../paper/images/overlaps_different_omega.pdf', bbox_inches='tight')
        plt.show()
    

    def plot_dipoles(self, time, dipole1, dipole2, dipole3, save_fig=False):
        plt.figure(figsize=(10,7))
        plt.plot(time* 0.5 * self.Omega / np.pi, dipole1, label=r'$\omega=8\Omega$')
        plt.plot(time* 0.5 * self.Omega / np.pi, dipole2, label=r'$\omega=4\Omega$')
        plt.plot(time* 0.5 * self.Omega / np.pi, dipole3, label=r'$\omega=32\Omega$')
        plt.xlabel(r'$t\Omega/2\pi}$', fontsize=22)
        plt.ylabel(r'$\overline{x}(t)$', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        plt.legend(fontsize=16)
        if save_fig==True:
            plt.savefig('../paper/images/dipoles_different_omega.pdf', bbox_inches='tight')
        plt.show()  
    
    def plot_fourier_analysis(self, time, dipole, overlap, xFFT, xfreqFFT, overlapFFT, overlapfreqFFT, save_dipole=False, save_overlap=False, save_fft_x=False, save_fft_overlap=False):    
        plt.figure(figsize=(16,5))
        plt.plot(time/np.max(time), dipole.real)
        plt.xlabel(r'$t/T_f$', fontsize=16)
        plt.ylabel(r'$\overline{x}(t)$', fontsize=16)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=18)
        plt.grid()
        if save_dipole==True:
            plt.savefig('../paper/images/dipole_laser_ONOFF.pdf', bbox_inches='tight')
        plt.show()
    
        plt.figure(figsize=(16,5))
        plt.plot(time/np.max(time), overlap.real)
        plt.xlabel(r'$t/T_f$', fontsize=16)
        plt.ylabel(r'$\xi (t)$', fontsize=16)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=18)
        plt.grid()
        if save_overlap==True:
            plt.savefig('../paper/images/overlap_laser_ONOFF.pdf', bbox_inches='tight')
        plt.show()

        plt.figure(figsize=(10,7))
        plt.plot(2*np.pi*xfreqFFT, np.abs(xFFT))
        plt.xlim([0, 3])
        plt.xlabel(r'$\omega_{fft}$', fontsize=22)
        plt.ylabel(r'$FFT[\overline{x}(t)]$', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        if save_fft_x==True:
            plt.savefig('../paper/images/x_spectrum.pdf', bbox_inches='tight')
        plt.show()
    
    
        plt.figure(figsize=(10,7))
        plt.plot(2*np.pi*overlapfreqFFT, np.abs(overlapFFT))
        plt.xlim([0, 3])
        plt.xlabel(r'$\omega_{fft}$', fontsize=22)
        plt.ylabel(r'$FFT[\xi_{T}(t)]$', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        if save_fft_overlap==True:
            plt.savefig('../paper/images/overlap_spectrum.pdf', bbox_inches='tight')
        plt.show()

    def plot_fft_dipole_comparison(self, freq1,freq2,freq3, fft1, fft2, fft3, save_fig=False):
        plt.figure(figsize=(10,7))
        freq1 = freq1[freq1>=0]
        freq2 = freq2[freq2>=0]
        freq3 = freq3[freq3>=0]
        fft1 = np.abs(fft1[0:len(freq1)])
        fft2 = np.abs(fft2[0:len(freq2)])
        fft3 = np.abs(fft3[0:len(freq3)])
        fft1 = fft1/np.max(fft1)
        fft2 = fft2/np.max(fft2)
        fft3 = fft3/np.max(fft3)

        
        plt.plot(2*np.pi*freq1, fft1, label=r'$\Omega=0.25, \omega=2.0$')
        plt.plot(2*np.pi*freq2, fft2, label=r'$\Omega=0.5, \omega=3.0$')
        plt.plot(2*np.pi*freq3, fft3, label=r'$\Omega=0.75, \omega=4.0$')
        plt.xlabel(r'$\omega_{fft}$', fontsize=22)
        plt.ylabel(r'$FFT[\overline{x}(t)]$', fontsize=22)
        plt.xlim([0, 2])
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        plt.legend(fontsize=14)

        if save_fig==True:
            plt.savefig('../paper/images/spectra_x_diff_omegas.pdf', bbox_inches='tight')
        plt.show()

    def plot_fft_overlap_comparison(self, freq1,freq2,freq3, fft1, fft2, fft3, save_fig=False):
        plt.figure(figsize=(10,7))
        freq1 = freq1[freq1>=0]
        freq2 = freq2[freq2>=0]
        freq3 = freq3[freq3>=0]
        fft1 = np.abs(fft1[0:len(freq1)])
        fft2 = np.abs(fft2[0:len(freq2)])
        fft3 = np.abs(fft3[0:len(freq3)])
        fft1 = fft1/np.max(fft1)
        fft2 = fft2/np.max(fft2)
        fft3 = fft3/np.max(fft3)
        
        plt.plot(2*np.pi*freq1, fft1, label=r'$\Omega=0.25, \omega=2.0$')
        plt.plot(2*np.pi*freq2, fft2, label=r'$\Omega=0.5, \omega=3.0$')
        plt.plot(2*np.pi*freq3, fft3, label=r'$\Omega=0.75, \omega=4.0$')
        plt.xlabel(r'$\omega_{fft}$', fontsize=22)
        plt.ylabel(r'$FFT[\xi_T(t)]$', fontsize=22)
        plt.xlim([0, 2])
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        plt.legend(fontsize=14)

        if save_fig==True:
            plt.savefig('../paper/images/spectra_xi_diff_omegas.pdf', bbox_inches='tight')
        plt.show()

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
    
    def eval_one_body_density(self, C, plot_ON=False, save_fig=False):
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
            if save_fig==True:
                plt.savefig('../paper/images/one_body_density_comparison.pdf', bbox_inches='tight')
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

        xFFT = np.fft.fft(dipole2)
        xfreqFFT = np.fft.fftfreq(len(dipole2), dt)
        overlapFFT = np.fft.fft(overlap2)
        overlapfreqFFT = np.fft.fftfreq(len(overlap2), dt)       
        
        return C2, time, dipole, overlap, xFFT, xfreqFFT, overlapFFT, overlapfreqFFT


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

    