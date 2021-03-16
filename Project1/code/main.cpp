#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <omp.h>

#include "System/system.h"
#include "Wavefunctions/gaussian.h"
#include "Wavefunctions/asymmetricGaussian.h"
#include "Hamiltonians/spherical.h"
#include "Hamiltonians/elliptical.h"
#include "Solvers/metropolis.h"
#include "Solvers/importanceSampling.h"
#include "Others/random_generator.h"
#include "Others/functions.h"
using namespace std;

int main(int argc, char *argv[]){

    // Select working mode : run "make" command. Then "./main x" where x is the number indicating working mode (see switch below)
    int selector = 0;
    if(argc>1){
        int a = stoi(argv[1]);
        assert(a>=0 && a<=7);
        selector = stoi(argv[1]);
    }

    // Information for the system
    const int dimension = 3;
    const int Nparticles = 3;

    // Information for the solvers
    const int Nsteps_final = (int) pow(2, 16); // MC steps for the final simulation
    const int NstepsThermal = (int) 1e5; // Fraction of septs to wait for the system thermalization
    const double step = 1.0; // only for metropolis
    const double D = 0.5; // only for importance sampling
    const double dt = 0.1; // only for importance sampling
    
    // Information for the hamiltonian
    const double a = 0; // Set the radius of the particles. a=0 is the non-interacting case
    double omegaXY = 1.0; // Only the elliptical hamiltonian distinguish between omegaXY and omegaZ
    double omegaZ = 1.0; 

    // Information for the wavefunction
    double alpha = 0.5; // variational parameter
    const double beta = 1.0; // Only for asymmetrical wavefunction
    
    // Others
    const double h = 1e-5; // Steplength for numerical derivatives and evaluations
    bool tofile = false; // Print on external file for resampling analysis (numerical methods

    // Parameters for the various type of simulations
    // Mode 1 - varying alpha
    const double alpha_min = 0.4; // in mode 1 (varying alpha) minimum alpha
    const double alpha_max = 0.6; // in mode 2 (varying alpha) maximum alpha
    const int N_alpha = 10; // in mode 1 (varying alpha) number of different alphas between alpha_min and alpha_max
    const bool alpha_to_file = false; // set true to save data to file

    // Mode 2 - varying dt
    const double dt_min = 1e-3; // in mode 2 (varying dt) minimum dt
    const double dt_max = 15; // in mode 2 (varying dt) maximum dt
    const int N_dt = 20; // in mode 2 (varying dt) number of different dts between dt_min dt_max
    const bool dt_to_file = false; // set true to save data to file

    // Mode 3 - varying N
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
    const double r_max = 3.0;
    int Nbins = 100;

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
    Functions functions(&system);

    // Choose options
    system.setHamiltonian(&spherical);
    system.setWavefunction(&gaussian);
    system.setSolver(&importance);
    system.setRandomGenerator(&randomgenerator);
    functions.printPresentation();
    
    // !!!!!!!! This should be more precise and should be put inside Function (time should be returned as another value)
    
    auto start = chrono::steady_clock::now(); // Store starting time to measure run time

    switch(selector){
        case 0: system.getSolver()->thermalize();
                functions.printResultsSolver(system.getSolver()->solve(false));  break; // Simple simulation
        case 1: system.getSolver()->thermalize();
                functions.printResultsSolver(system.getSolver()->solve(h)); break; // Simple simulation with numerical derivative
        case 2: functions.solve_varying_alpha(alpha_min, alpha_max, N_alpha, alpha_to_file); break;
        case 3: functions.solve_varying_dt(dt_min, dt_max, N_dt, dt_to_file); break;
        case 4: functions.solve_varying_N(Ns, N_to_file); break;
        case 5: best_alpha = functions.gradientDescent(initial_alpha, gamma, tolerance, Nmax_gradient, Nsteps_gradient); 
                cout << scientific << setprecision(5) << "best alpha= " << best_alpha << endl; break;
        case 6: system.getSolver()->thermalize();
                functions.printResultsSolver(system.getSolver()->solve(r_max, Nbins)); break;
        case 7: 
                int Nthreads = (int) omp_get_max_threads();
                int Ni = (int) Nsteps_final/Nthreads;
                #pragma omp parallel for schedule(static) num_threads(Nthreads) \
                shared(Ni, Nthreads)
                for(int i=0; i<Nthreads; i++){

                        System sys(dimension, Nparticles);

                        Elliptical ellipt(&sys, omegaXY, omegaZ);
                        AsymmetricGaussian asymmgauss(&sys, alpha, beta, a);
                        ImportanceSampling imp(&sys, Nsteps_final, NstepsThermal, dt, D, tofile);
                        RandomGenerator randomgen;
                        Functions funcs(&sys);

                        sys.setHamiltonian(&ellipt);
                        sys.setWavefunction(&asymmgauss);
                        sys.setSolver(&imp);
                        sys.setRandomGenerator(&randomgen);
                        sys.getSolver()->thermalize();

                        cout << sys.getSolver()->solve(false)[0] << endl;  
                }
                break;
    }
    
    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << fixed << "Simulation termined. Total running time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl; // Print run time*/
}