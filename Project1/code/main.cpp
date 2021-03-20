#include "Others/functions.h"
#include "System/system.h"
#include "Wavefunctions/gaussian.h"
#include "Wavefunctions/asymmetricGaussian.h"
#include "Hamiltonians/spherical.h"
#include "Hamiltonians/elliptical.h"
#include "Solvers/metropolis.h"
#include "Solvers/importanceSampling.h"
#include "Others/random_generator.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]){

    // Select working mode : run "make" command. Then "./main x" where x is the number indicating working mode (see switch below)
    int selector = 0;
    bool parallel = false;
    if(argc>1){
        parallel = (bool) stoi(argv[1]);
        int a = stoi(argv[2]);
        assert(a>=0 && a<=6);
        selector = stoi(argv[2]);
    }

    // Information for the system
    const int dimension = 3;
    const int Nparticles = 10;

    // Information for the solvers
    const int Nsteps_final = (int) pow(2,21); // MC steps for the final simulation
    const int NstepsThermal = (int) 1e5; // Fraction of septs to wait for the system thermalization
    const double step = 1.0; // only for metropolis
    const double D = 0.5; // only for importance sampling
    const double dt = 0.01; // only for importance sampling
    
    // Information for the hamiltonian
    const double a = 0.0043; // Set the radius of the particles. a=0 is the non-interacting case
    double omegaXY = 1.0; // Only the elliptical hamiltonian distinguish between omegaXY and omegaZ
    double omegaZ = 2.82843; 

    // Information for the wavefunction
    double alpha = 0.5; // variational parameter
    const double beta = 2.82843; // Only for asymmetrical wavefunction
    
    // Others
    bool tofile = true; // Print on external file for resampling analysis (numerical methods)

    // Parameters for the various type of simulations

    // Mode 0 - normal

    // Mode 1 - numerical
    const double h = 1e-5; // Steplength for numerical derivatives and evaluations

    // Mode 2 - varying alpha
    const double alpha_min = 0.4; // in mode 1 (varying alpha) minimum alpha
    const double alpha_max = 0.6; // in mode 2 (varying alpha) maximum alpha
    const int N_alpha = 10; // in mode 1 (varying alpha) number of different alphas between alpha_min and alpha_max
    const bool alpha_to_file = false; // set true to save data to file

    // Mode 3 - varying dt
    const double dt_min = 1e-3; // in mode 2 (varying dt) minimum dt
    const double dt_max = 10; // in mode 2 (varying dt) maximum dt
    const int N_dt = 10; // in mode 2 (varying dt) number of different dts between dt_min dt_max
    const bool dt_to_file = false; // set true to save data to file

    // Mode 4 - varying N
    vector<int> Ns {5, 10, 15}; // in mode 3 (varying N) different values of N
    const bool N_to_file = false; // set true to save data to file

    // Mode 5 - Gradient Descent
    double best_alpha = 0.0;
    double initial_alpha = 0.49;
    const double gamma = 2e-3; // Learning rate
    const unsigned int Nmax_gradient = 100; // Max steps alowed
    const unsigned int Nsteps_gradient = (int) 1e5; // MC steps during the gradient descent
    const double tolerance = 1e-8; // Conditoin to stop the gradient descent

    // Mode 6 - One Body density
    const double r_max = 4.0;
    const int Nbins = 200;




    // BOOL PER SELEZIONI VARIE
    bool use_elliptic = true;
    bool use_asymmetric = true;
    bool use_importance = true;    
        
    // !!!!!!!! This should be more precise and should be put inside Function (time should be returned as another value)
    
    auto start = chrono::steady_clock::now(); // Store starting time to measure run time

    if(parallel){
        int Nthreads = omp_get_max_threads();
        int Ni = (int) Nsteps_final/Nthreads;
        //!!!!!!! omp_get_thread_num();
        #pragma omp parallel for num_threads(Nthreads) schedule(static, 1) shared(Ni, Nthreads)
        for(int i=0; i<Nthreads; i++){
            System system(dimension, Nparticles);

            // Hamiltonians
            Spherical spherical(&system, omegaXY);
            Elliptical elliptical(&system, omegaXY, omegaZ);

            // Wavefunctions
            Gaussian gaussian(&system, alpha);
            AsymmetricGaussian asymmgaussian(&system, alpha, beta, a);

            // Solvers
            Metropolis metropolis(&system, Ni, NstepsThermal, step, tofile);
            ImportanceSampling importance(&system, Ni, NstepsThermal, dt, D, tofile);
        
            // Others
            RandomGenerator randomgenerator;
            Functions functions(&system, true);

            // Choose options
            if(use_elliptic) system.setHamiltonian(&elliptical);
            else system.setHamiltonian(&spherical);

            if(use_asymmetric) system.setWavefunction(&asymmgaussian);
            else system.setWavefunction(&gaussian);

            if(use_importance) system.setSolver(&importance);
            else system.setSolver(&metropolis);
            
            system.setRandomGenerator(&randomgenerator);
            if(omp_get_thread_num()==0) {functions.printPresentation();}

            switch(selector){
                case 0: functions.solve_singleRun();  break; // Simple simulation
                case 1: functions.solve_singleRun(h); break; // Simple simulation with numerical derivative
                case 2: functions.solve_varying_alpha(alpha_min, alpha_max, N_alpha, alpha_to_file); break;
                case 3: functions.solve_varying_dt(dt_min, dt_max, N_dt, dt_to_file); break;
                case 4: functions.solve_varying_N(Ns, N_to_file); break;
                case 5: functions.gradientDescent(initial_alpha, gamma, tolerance, Nmax_gradient, Nsteps_gradient); break;
                case 6: functions.solve_singleRun(r_max, Nbins); break;
            }

        }
    
    }   else    {
        
        System system(dimension, Nparticles);

        // Hamiltonians
        Spherical spherical(&system, omegaXY);
        Elliptical elliptical(&system, omegaXY, omegaZ);

        // Wavefunctions
        Gaussian gaussian(&system, alpha);
        AsymmetricGaussian asymmgaussian(&system, alpha, beta, a);

        // Solvers
        Metropolis metropolis(&system, Nsteps_final, NstepsThermal, step, tofile);
        ImportanceSampling importance(&system, Nsteps_final, NstepsThermal, dt, D, tofile);

        // Others
        RandomGenerator randomgenerator;
        Functions functions(&system, false);

        // Choose options
        if(use_elliptic) system.setHamiltonian(&elliptical);
        else system.setHamiltonian(&spherical);

        if(use_asymmetric) system.setWavefunction(&asymmgaussian);
        else system.setWavefunction(&gaussian);

        if(use_importance) system.setSolver(&importance);
        else system.setSolver(&metropolis);
            
        system.setRandomGenerator(&randomgenerator);
        functions.printPresentation();

        switch(selector){
        case 0: functions.solve_singleRun();  break; // Simple simulation
        case 1: functions.solve_singleRun(h); break; // Simple simulation with numerical derivative
        case 2: functions.solve_varying_alpha(alpha_min, alpha_max, N_alpha, alpha_to_file); break;
        case 3: functions.solve_varying_dt(dt_min, dt_max, N_dt, dt_to_file); break;
        case 4: functions.solve_varying_N(Ns, N_to_file); break;
        case 5: functions.gradientDescent(initial_alpha, gamma, tolerance, Nmax_gradient, Nsteps_gradient); break;
        case 6: functions.solve_singleRun(r_max, Nbins); break;
        }


    }

    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << fixed << "Simulation termined. Total running time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl; // Print run time*/
    
}

