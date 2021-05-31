from logging import log
from operator import pos
from os import TMP_MAX, system
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import shape
from numpy.lib.type_check import imag
import scipy.linalg, scipy.integrate
from quantum_systems import ODQD, GeneralOrbitalSystem
from IPython.display import display, clear_output
import time as timeimp
from matplotlib.animation import FuncAnimation, PillowWriter
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

class GHF:
    """
        Initializes a system in the general spin representation. The default values in the creator are referred to the main analyzed case.
        A more detailed description is provided in the paper folder at https://github.com/Matteo294/FYS4411.
        The quantum-system documentation appears in https://schoyen.github.io/quantum-systems/

        **Args**:
        
            l (int) : number of harmonic oscillator eigenstates used for the expansion of the radial part of the single-particle wavefunctions

            grid_length (int) : extension of the mesh goes from -grid_length to grid_length

            num_grid_points (int) : number of points in the mesh

            alpha (float) : strength parameter in the shielded Coulomb potential

            a (float) : shielding parameter in the shielded Coulomb potential

            Omega (float) : frequency of the harmonic oscillator potential in which the electrons are trapped

            omega (float) : frequency of the laser source

            epsilon0 (float) : amplitude of the sinusoidal potential associated to the laser

            nparticles (int) : number of electrons in the system

            antisymmetrize (bool) : if True, antisymmetrizes the system.u matrix

        **Attributes**:
        
            potential (func) : ODQD.HOPotential(Omega) 

            system : GeneralOrbitalSystem object (see https://schoyen.github.io/quantum-systems/)

            Omega (float) : frequency of the harmonic oscillator potential in which the electrons are trapped

            omega (float) : frequency of the laser source

            epsilon0 (float) : amplitude of the sinusoidal potential associated to the laser

            nparticles (int) : number of electrons in the system

    """
    def __init__(self, l=10, grid_length=10, num_grid_points=201, alpha=1.0, a=0.25, Omega=0.25, omega=2, epsilon0=1, nparticles=2, antisymmetrize=True):
        
        self.potential = ODQD.HOPotential(Omega)
        
        odho = ODQD(l, grid_length, num_grid_points, a, alpha, potential=self.potential)
        self.system = GeneralOrbitalSystem(n=nparticles, basis_set=odho, anti_symmetrize=antisymmetrize)
        
        self.nparticles = nparticles
        self.Omega = Omega
        self.omega = omega
        self.epsilon0 = epsilon0
        
    
    def laser_potential(self, t):
        """
        Returns the value of the laser potential at time t.

        **Args**:
            t (float) : time instant
        
        **Returns**:
            epsilon0 * sin( omega * t ) (float) 

        """
        return self.epsilon0 * np.sin( self.omega * t )
    
    def fill_density_matrix(self, C):
        r"""
        Returns the density matrix evaluated using the coefficient matrix C, according to
        
        :math:`M_{\alpha\beta} = \sum_i^{occ} C_{\alpha,i}^* C_{\beta,i}`

        **Args**:
            C (np.ndarray) : coefficient matrix
        
        **Returns**:
            density_matrix (np.ndarray)

        """
        density = np.zeros( C.shape, dtype=np.complex128)
        for i in range(self.nparticles):
            density += np.outer( np.conjugate(C[:,i]), C[:,i])
        
        return density
        
    def fill_fock_matrix(self, C, t):
        r"""
        Returns the Fock matrix evaluated using the coefficient matrix C evaluated at time t, according to

        :math:`f_{\mu\nu}(t) =  h_{\mu\nu}^{ho} + x_{\mu\nu} \varepsilon_0 \sin(\omega t) + \sum_{j}^{occ} \sum_{\gamma\delta} C_{\gamma,j}^* C_{\delta,j} u^{\mu\gamma}_{\nu\delta,AS}`

        **Args**:
            C (np.ndarray) : coefficient matrix

            t (float) : time instant
        
        **Returns**:
            fock_matrix (np.ndarray)

        """
        f = np.zeros(self.system.h.shape, dtype=np.complex128)    
        density = self.fill_density_matrix(C)
        
        f = np.einsum('ij,aibj->ab', density, self.system.u, dtype=np.complex128)
        f += self.system.h
        f += self.system.position[0] * self.laser_potential(t)
        return f
    
    def evaluate_total_energy(self, C):
        r"""
        Returns the total energy of the system, given by
        
        :math:`E[H] = \sum_{i=1}^2 \sum_{\alpha,\beta} C_{\alpha,i}^* C_{\beta,i} h_{\alpha\beta}^{ho} + \frac{1}{2} \sum_{i,j=1}^2 \sum_{\alpha\beta\gamma\delta} C_{\alpha,i}^* C_{\gamma,j}^* C_{\beta,i} C_{\delta,j}  u^{\alpha\gamma}_{\beta\delta,AS}`

        **Args**:
            C (np.ndarray) : coefficient matrix
        
        **Returns**:
            energy (complex)

        """
        density = self.fill_density_matrix(C)
        rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
        energy = np.einsum('ij,ij', density, self.system.h) + 0.5*np.einsum('abcd,acbd', rhorho, self.system.u, dtype=np.complex128)
        
        return energy
    
    def eval_one_body_density(self, C):
        r"""
        Returns the one-body density of the system, according to

        :math:`\rho(x) = \sum_i^{occ} \sum_{\alpha\beta} C_{\alpha,i}^* C_{\beta,i} \chi_{\alpha}^*(x) \chi_{\beta}(x)`

        **Args**:
            C (np.ndarray) : coefficient matrix

        **Returns**:
            one_body_density (np.array)

        """
        obd = np.zeros( len(self.system.grid) )
        density = self.fill_density_matrix(C)
        obd = np.einsum('mi,mn,ni->i', self.system.spf, density, self.system.spf, dtype=np.complex128)

        return obd
    
    def eval_dipole(self, C):
        r"""
        Returns the expected value of the position operator for the system, given by

        :math:`\overline{x}(t) = \sum_i^{occ} \sum_{\alpha\beta} C_{\alpha,i}^*(t) C_{\beta,i} (t) \langle \chi_{\alpha}(x) \vert \hat{x} \vert \chi_{\beta}(x)\rangle`

        **Args**:
            C (np.ndarray) : coefficient matrix
        
        **Returns**:            
            dipole (complex) : expected value for position operator

        """
        density = self.fill_density_matrix(C)
        dipole = np.einsum('mn,mn', density, self.system.position[0])
        return dipole
    
    def solve_TIHF(self, tolerance, max_iter, print_ON=False, eval_energy_per_step=False, eval_delta=False):
        r"""
        Solves the Ruthaan-Hall equations for the system.

        **Args**:
            tolerance (float) : stopping condition for the algorithm.
            
            max_iter (int) : maximum number of iterations before the algorithm stops

            print_ON (bool) : if True prints if the convergence or the max number of iterations has been reached.

            energy_per_step (bool) : if True, the energy of the system is evaluated at every step

            eval_delta (bool) : if True, the :math:`\Delta` parameter of the system is evaluated at every step
        
        **Returns**:
            epsilon (np.array) : final eigenvalues
            
            C (np.ndarray) : final coefficient matrix
            
            energy_per_step (np.array) : None if energy_per_step==False
            
            delta_per_step (np.array) : None if delta_per_step==False

        """
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
        r"""
            Returns the right-hand side of the Ruthaan-Hall equations with the laser source on, corresponding to

            :math:`\dot{C}(t) = -i f(t)C(t)`
            
            The output is reshaped into an array in order to be used into self.solve_TDHF()

            **Args**:
                t (float) : time

                C (np.ndarray) : coefficient matrix
            
            **Returns**:
                rhs (np.ndarray)
        
        """
        C = np.reshape(C, self.system.h.shape)
        f = self.fill_fock_matrix(C, t)
        res = -1.0j * (f @ C)
        res = np.reshape(res, (len(res)**2) )
        return res
    
    def rhsf_OFF(self, t, C):
        r"""
            Returns the right-hand side of the Ruthaan-Hall equations with the laser source off, given by

            :math:`\dot{C}(t) = -i f(t)C(t)`

            The output is reshaped into an array in order to be used into self.solve_TDHF()

            **Args**:
                t (float) : time

                C (np.ndarray) : coefficient matrix
            
            **Returns**:
                rhs (np.ndarray)

        """
        C = np.reshape(C, self.system.h.shape)
        f = self.fill_fock_matrix(C, t=0)
        res = -1.0j * (f @ C)
        res = np.reshape(res, (len(res)**2) )
        return res
    
    def solve_TDHF(self, tstart, dt, t_max, C0, eval_overlap=False, eval_dipole=False, eval_energy=False, laser_ON=True):
        r"""
        Solves iteratively the time-dependent Ruthaan-Hall equations for the system.

        **Args**:
            t_start (float) : time at which the evolution begins

            dt (float) : time step

            t_max (float) : time at which the evolution ends

            C0 (np.ndarray) : coefficient matrix at time t_start

            eval_overlap (bool) : if True, the overlap with the wavefunction at t=t_start is evaluated at every instant

            eval_dipole (bool) : if True, the expected value for the dipole operator is evaluated at every instant
            
            eval_energy (bool) : if True, the total energy is evaluated at every instant

            laser_ON (bool) : if True, the laser source is on
        
        **Return**: 
            C (np.ndarray) : coefficient matrix at t=t_max
            
            time (np.array) : time instants from tstart to t_max spaced by dt
            
            overlap (np.array) : None if eval_overlap==False
            
            dipole (np.array) : None if eval_dipole==False
            
            energy (np.array) : None if eval_energy==False

        """
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
        if eval_energy==True:
            energy = np.zeros( (i_max), dtype=np.complex128 )
        else:
            energy=None

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
            if eval_energy==True:
                energy[i] = self.evaluate_total_energy(C)
            i += 1
            if i%10==0:
                print('\rprogress = %.2f' % (100*i/i_max), end='')
        
        print('\r', end='')
        C = np.array(C)

        return C, time, overlap, dipole, energy

    def fourier_analysis(self, tolerance, max_iter, t_laser_ON, t_max, dt, eval_energy=False):
        r"""
        Solves the time-independent Ruthaan-Hall equations and then performs a time-evolution of the system switching off the laser source at a certain time.
        Performs the Fourier analysis on the curves for overlap and dipole moment obtained for t>t_laser_ON

        **Args**:
            tolerance (float) : stopping condition for the time-independent solver.
            
            max_iter (int) : maximum number of iterations before the time-independent solver stops.

            t_laser_ON (float) : instant at which the laser is switched off

            t_max (float) : final instant for time evolution

            dt (float) : time step

            eval_energy (bool) : leave false for faster computation

        **Returns**:
            C2 (np.ndarray) : coefficient matrix for t=t_max

            time (np.array) : time instants from tstart to t_max spaced by dt

            dipole (np.array) : dipole evaluated at each time instant

            overlap (np.array) : overlap evaluated at each time instant

            dipoleFFT (np.array) : fft of dipole values for t>t_laser_ON

            dipolefreqFFT (np.array) : array of frequencies corresponding to dipoleFFT

            overlapFFT (np.array) : fft of overlap values for t>t_laser_ON

            overlapfreqFFT (np.array) : array of frequencies corresponding to overlapFFT

            energy (np.array) : total energy at each time instant

        """
        epsilon, C0, energy_per_step, delta_per_step = self.solve_TIHF(tolerance=tolerance, max_iter=max_iter, print_ON=False)
        C1, time1, overlap1, dipole1, energy1 = self.solve_TDHF(0, dt, t_laser_ON, C0, eval_overlap=True, eval_dipole=True, eval_energy=eval_energy, laser_ON=True)
        C2, time2, overlap2, dipole2, energy2 = self.solve_TDHF(t_laser_ON,  dt, t_max, C1, eval_overlap=True, eval_dipole=True, eval_energy=eval_energy, laser_ON=False)
        time = np.concatenate((time1,time2))
        dipole = np.concatenate((dipole1, dipole2))
        overlap = np.concatenate((overlap1, overlap2))
        if eval_energy==True:
            energy = np.concatenate((energy1, energy2))
        else:
            energy = None

        xFFT = np.fft.fft(dipole2)
        xfreqFFT = np.fft.fftfreq(len(dipole2), dt)
        overlapFFT = np.fft.fft(overlap2)
        overlapfreqFFT = np.fft.fftfreq(len(overlap2), dt)       
        
        return C2, time, dipole, overlap, xFFT, xfreqFFT, overlapFFT, overlapfreqFFT, energy


    def gif_generator(self, dt, t_max, C0, save_every_n=500):
        r"""
        Generates an animated image (.gif) with the time evolution of the one-body density and the time evolution of the one body potential.

        **Args**:
            dt (float) : time step

            t_max (float) : total duration of the time evolution

            C0 (np.ndarray) : coefficient matrix at time t=0

            save_every_n (int) : number of time steps between two successive frames in the final .gif file

        """
        fig, ax = plt.subplots(figsize=(10,7))
        ax.set_xlim([-10, 10])
        ax.set_ylim([-0.1, 1.5])
        ax.set_xlabel(r'$x$', fontsize=22)
        ax.tick_params(axis='both', which='major', pad=5, labelsize=18)
        ax.grid()

        epsilon, C0, energy_per_step, delta_per_step = self.solve_TIHF(tolerance=1e-6, max_iter=100)
        line1, = ax.plot(self.system.grid, self.eval_one_body_density(C0).real)
        line2, = ax.plot(self.system.grid, self.potential(self.system.grid))
        props = dict(boxstyle='round', facecolor='white')
        text = ax.text(-9, 1.4, r'$t/T_{max}=0$', fontsize=16, bbox=props, animated=True)

        integrator = scipy.integrate.ode(self.rhsf).set_integrator('zvode')
        integrator.set_initial_value( np.reshape(C0, len(C0)**2 ), 0)
        
        i_max = i_max = int( np.floor(t_max / dt / save_every_n) )

        anim = FuncAnimation(fig, self.anima, fargs=(text, line1, line2, integrator, dt, t_max, save_every_n), frames=i_max)
        writer = PillowWriter(fps=20)  
        anim.save("animation.gif", writer=writer)  
    

    def anima(self, i, text, line1, line2, integrator, dt, t_max, save_every_n):
        r"""
        Iterative function needed for self.gif_generator

        """
        for i in range(save_every_n):
            C = integrator.integrate(integrator.t+dt)
            C = np.matrix(np.reshape(C, self.system.h.shape))

        line1.set_ydata(self.eval_one_body_density(C).real)
        line2.set_ydata(self.potential(self.system.grid) + self.system.grid*self.laser_potential(integrator.t))
        text.set_text(r'$t/T_{max}=$'+ f'{integrator.t/t_max:.2}')
        
    
    def plot_AO(self): 
        r"""
            Plots the square modulus of the Atomic Orbitals (harmonic oscillator eigenstates).
            The zero-level for each single particle state is its corresponding eigenvalue.

        """
        plt.figure(figsize=(10,7))
        plt.plot(self.system.grid, self.potential(self.system.grid))

        for i in range(int(self.system.l/2)):
            plt.plot( self.system.grid, np.abs(self.system.spf[2*i]) ** 2 + self.system.h[2*i, 2*i].real, label=r"$\chi_{" + f"{i}" + r"}$" )

        plt.grid()
        plt.legend()
        plt.show()
    
    def plot_MO(self, epsilon, C):
        r"""
            Plots the square modulus of the Molecular Orbitals (after the basis change).
            The zero-level for each single particle state is its corresponding eigenvalue.

        """
        plt.figure(figsize=(10,7))
        plt.plot(self.system.grid, self.potential(self.system.grid))
        
        for i in range(self.system.l):
            to_plot = np.zeros( self.system.grid.shape, np.complex128 )
            for j in range(self.system.l):
                to_plot += self.system.spf[j] * C[j,i]
            plt.plot( self.system.grid,  np.abs(to_plot)**2 + epsilon[i], label=r"$\phi_{" + f"{i}" + r"}$" )

        plt.grid()
        plt.legend()
        plt.show()
    
    def plot_one_body_density(self, one_body_density, save_fig=False):
        r"""
            Plots the comparison between the one-body density and the result by Zanghellini et al..

            **Args**:
                one_body_density (np.array) : result of self.eval_one_body_density()

        """
        plt.figure(figsize=(13, 13))
        img = plt.imread("theoretical_density.png")
        ext = [-6.0, 6.0, 0.00, 0.4]
        plt.imshow(img, zorder=0, extent=ext)
        aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
        plt.gca().set_aspect(aspect)
        plt.plot(self.system.grid , one_body_density.real)
        plt.xlabel('x', fontsize=16)
        plt.ylabel(r'$\rho(x)$', fontsize=16)
        plt.xlim([-10, 10])
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
        plt.grid()
        if save_fig==True:
            plt.savefig('../paper/images/one_body_density_comparison.pdf', bbox_inches='tight')
        plt.show()
    
    def plot_overlap(self, time, overlap, save_fig=False):
        r"""
            Plots the comparison between the overlaps evaluated in the time-dependent solver and by Zanghellini et al..

            **Args**:
                time (np.array) : array of time instants

                overlap (np.array) : overlap for each time instant
        """
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
    
    def plot_fourier_analysis(self, time, dipole, overlap, xFFT, xfreqFFT, overlapFFT, overlapfreqFFT, energy=None, save_dipole=False, save_overlap=False, save_energy=False, save_fft_x=False, save_fft_overlap=False):    
        r"""
            Plots the results of the Fourier analysis.

            **Args**:
                time (np.array) : array of time instants

                dipole (np.array) : dipole for each time instant

                overlap (np.array) : overlap for each time instant
                
                xFFT (np.array) : fast Fourier transform of the dipole signal
                
                xfreqFFT (np.array) : frequency spectrum of the dipole signal (useful for plotting) 

                overlapFFT (np.array) : fast Fourier transform of the overlap signal
                
                overlapfreqFFT (np.array) : frequency spectrum of the overlap signal (useful for plotting) 

                energy (np.array) : energy for each time instant

        """
        plt.figure(figsize=(16,5))
        plt.plot(time/np.max(time), dipole.real)
        plt.xlabel(r'$t/T_f$', fontsize=16)
        plt.ylabel(r'$\overline{x}(t)$ [a.u]', fontsize=16)
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

        if energy is not None:
            plt.figure(figsize=(16,5))
            plt.plot(time/np.max(time), energy.real)
            plt.xlabel(r'$t/T_f$', fontsize=16)
            plt.ylabel(r'$Energy (t)$', fontsize=16)
            ax = plt.gca()
            ax.tick_params(axis='both', which='major', pad=5, labelsize=18)
            plt.grid()
            if save_energy==True:
                plt.savefig('../paper/images/energy_laser_ONOFF.pdf', bbox_inches='tight')
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
    
    
    def plot_onebody_comparison(self, onebodyg, onebodyr, save_fig=False):
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
        plt.figure(figsize=(16, 16))
        img = plt.imread("theoretical_density.png")
        ext = [-6.0, 6.0, 0.00, 0.4]
        plt.imshow(img, zorder=0, extent=ext)
        aspect = img.shape[0]/float(img.shape[1])*((ext[1]-ext[0])/(ext[3]-ext[2]))
        plt.gca().set_aspect(aspect)
        plt.plot(self.system.grid, onebodyg.real, label='GHF')
        plt.plot(self.system.grid, onebodyr.real, label='RHF')
        plt.xlabel(r'$x \; [a.u.]$', fontsize=16)
        plt.ylabel(r'$\rho(x)$', fontsize=16)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=16)
        plt.grid()
        plt.legend(fontsize=12)
        if save_fig==True:
            plt.savefig('../paper/images/onebody_density_comp_article.pdf', bbox_inches='tight')
        plt.show()
    
    def plot_energy_per_iteration(self, energy_ghf, energy_rhf, save_fig=False):
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
        plt.figure(figsize=(10,7))
        plt.plot(np.arange(len(energy_ghf)), energy_ghf.real, label='GHF')
        plt.plot(np.arange(len(energy_rhf)), energy_rhf.real, label='RHF')
        plt.xlabel('Iteration', fontsize=22)
        plt.ylabel('Energy [Hartree]', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        plt.legend(fontsize=14)

        if save_fig==True:
            plt.savefig('../paper/images/energy_at_every_step.pdf', bbox_inches='tight')
        plt.show()
        
    def plot_delta_per_iteration(self, delta_ghf, delta_rhf, tolerance, save_fig=False):
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
        plt.figure(figsize=(10,7))
        plt.plot(np.arange(len(delta_ghf)), delta_ghf.real, label='GHF')
        plt.plot(np.arange(len(delta_rhf)), delta_rhf.real, label='RHF')
        plt.plot(np.arange(len(delta_rhf)), tolerance*np.ones(len(delta_ghf)), linestyle='--', color='red', label=r'$\delta$')
        plt.xlabel('Iteration', fontsize=22)
        plt.ylabel(r'$\Delta$ [Hartree]', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        ax.set_yscale('log')
        plt.grid()
        plt.legend(fontsize=14)

        if save_fig==True:
            plt.savefig('../paper/images/delta_at_every_step.pdf', bbox_inches='tight')
        plt.show()

    def plot_dipole(self, time, dipole, save_fig=False):
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
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
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
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
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
        plt.figure(figsize=(10,7))
        plt.plot(time* 0.5 * self.Omega / np.pi, dipole1, label=r'$\omega=8\Omega$')
        plt.plot(time* 0.5 * self.Omega / np.pi, dipole2, label=r'$\omega=4\Omega$')
        plt.plot(time* 0.5 * self.Omega / np.pi, dipole3, label=r'$\omega=32\Omega$')
        plt.xlabel(r'$t\Omega/2\pi}$', fontsize=22)
        plt.ylabel(r'$\overline{x}(t)$ [a.u.]', fontsize=22)
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', pad=5, labelsize=22)
        plt.grid()
        plt.legend(fontsize=16)
        if save_fig==True:
            plt.savefig('../paper/images/dipoles_different_omega.pdf', bbox_inches='tight')
        plt.show()  

    def plot_fft_dipole_comparison(self, freq1,freq2,freq3, fft1, fft2, fft3, save_fig=False):
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
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
        r"""
            FACCIAMO SENZA DOCUMENTAZIONE DAI
        """
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
        
        
    
    



class RHF(GHF):
    r"""
        Initializes a system in the restricted spin representation. The default values in the creator are referred to the analyzed case.
        A more detailed description is provided in the paper folder in this same GitHub repository.
        The quantum-system documentation appears in https://schoyen.github.io/quantum-systems/

        **Args**:
            l (int) : number of harmonic oscillator eigenstates used for the expansion of the radial part of the single-particle wavefunctions

            grid_length (int) : extension of the mesh goes from -grid_length to grid_length

            num_grid_points (int) : number of points in the mesh

            alpha (float) : strength parameter in the shielded Coulomb potential

            a (float) : shielding parameter in the shielded Coulomb potential

            Omega (float) : frequency of the harmonic oscillator potential in which the electrons are trapped

            omega (float) : frequency of the laser source

            epsilon0 (float) : amplitude of the sinusoidal potential associated to the laser

            nparticles (int) : number of electrons in the system

            antisimmetrize (bool) : MUST remain False

        **Attributes**:
            potential (func) : ODQD.HOPotential(Omega)  

            system : ODQD object (see https://schoyen.github.io/quantum-systems/)

            Omega (float) : frequency of the harmonic oscillator potential in which the electrons are trapped

            omega (float) : frequency of the laser source

            epsilon0 (float) : amplitude of the sinusoidal potential associated to the laser

            nparticles (int) : number of electrons in the system

    """

    def __init__(self, l=10, grid_length=10, num_grid_points=201, alpha=1.0, a=0.25, Omega=0.25, omega=2, epsilon0=1, nparticles=2, potential=None, antisymmetrize=False):
        
        self.potential = ODQD.HOPotential(Omega)
                
        self.system = ODQD(l, grid_length, num_grid_points, a, alpha, potential=self.potential)
        self.nparticles = nparticles
        self.Omega = Omega
        self.omega = omega
        self.epsilon0 = epsilon0
    
    def fill_density_matrix(self, C):
        r"""
        Returns the density matrix evaluated using the coefficient matrix C, according to
        
        :math:`M_{\alpha\beta} = \sum_i^{nparticles/2} C_{\alpha,i}^* C_{\beta,i}`

        **Args**:
            C (np.ndarray) : coefficient matrix
        
        **Returns**:
            density_matrix (np.ndarray)

        """

        density = np.zeros( C.shape, dtype=np.complex128)
        for i in range(int(self.nparticles/2)):
            density += np.outer( np.conjugate(C[:,i]), C[:,i])
        
        return density
        return 

    def fill_fock_matrix(self, C, t):
        r"""
        Returns the Fock matrix evaluated using the coefficient matrix C evaluated at time t, according to

        :math:`f_{\mu\nu} =  h_{\mu\nu}^{ho} + 2\sum_{j}^{n/2} \sum_{\gamma\delta} C_{\gamma,j}^* C_{\delta,j} u^{\mu\gamma}_{\nu\delta} - \sum_{j}^{n/2} \sum_{\gamma\delta} C_{\gamma,j}^* C_{\delta,j} u^{\mu\gamma}_{\delta\nu}`

        **Args**:
            C (np.ndarray) : coefficient matrix

            t (float) : time instant
        
        **Returns**:
            fock_matrix (np.ndarray)

        """
        f = np.zeros(self.system.h.shape, dtype=np.complex128)    
        density = self.fill_density_matrix(C)
        f = 2*np.einsum('ij,aibj->ab', density, self.system.u, dtype=np.complex128)
        f -= np.einsum('ij,aijb->ab', density, self.system.u, dtype=np.complex128)
        f += self.system.h
        f += self.system.position[0] * self.laser_potential(t)

        return f
    
    def evaluate_total_energy(self, C):
        r"""
        Returns the total energy of the system, given by
        
        :math:`E[H] = 2\sum_{i=1}^{n/2} \sum_{\alpha,\beta} C_{\alpha,i}^* C_{\beta,i} h_{\alpha\beta}^{ho} + 2 \sum_{i,j=1}^{n/2} \sum_{\alpha\beta\gamma\delta} C_{\alpha,i}^* C_{\gamma,j}^* C_{\beta,i} C_{\delta,j}  u^{\alpha\gamma}_{\beta\delta} - \sum_{i,j=1}^{n/2} \sum_{\alpha\beta\gamma\delta} C_{\alpha,i}^* C_{\gamma,j}^* C_{\beta,i} C_{\delta,j}  u^{\alpha\gamma}_{\delta\beta}`

        **Args**:
            C (np.ndarray) : coefficient matrix
        
        **Returns**:
            energy (complex)

        """
        density = self.fill_density_matrix(C)
        rhorho = np.einsum('ab,cd->abcd', density, density, dtype=np.complex128)
        energy = 2*np.einsum('ij,ij', density, self.system.h) 
        energy += 2*np.einsum('abcd,acbd', rhorho, self.system.u, dtype=np.complex128)
        energy -= np.einsum('abcd,acdb', rhorho, self.system.u, dtype=np.complex128)
        
        return energy
    
    def eval_one_body_density(self, C):
        r"""
        Returns the one-body density of the system, according to

        :math:`\rho(x) = 2 \sum_i^{n/2} \sum_{\alpha\beta} C_{\alpha,i}^* C_{\beta,i} \chi_{\alpha}^*(x) \chi_{\beta}(x)`

        **Args**:
            C (np.ndarray) : coefficient matrix

        **Returns**:
            one_body_density (np.array)

        """
        obd = np.zeros( len(self.system.grid) )
        density = self.fill_density_matrix(C)
        obd = 2*np.einsum('mi,mn,ni->i', self.system.spf, density, self.system.spf, dtype=np.complex128)

        return obd

    
