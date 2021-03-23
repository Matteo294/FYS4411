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


// Set to 1 to use the parallelization
#define RUN_PARALLEL 1
// Set the correct number depending on your number of cores
#define NTHREADS 4

#define DIMENSION 3
#define NPARTICLES 50
#define USE_ASYMMETRIC 1
#define USE_ELLIPTICAL 1
#define USE_IMPORTANCE 1
#define TO_FILE 0

using namespace std;

int main(int argc, char *argv[]){

    // adjustable parameters
    const int Nsteps_final = (int) pow(2,21); // MC steps for the final simulation
    const int NstepsThermal = (int) 1e5; // Fraction of septs to wait for the system thermalization
    double alpha = 0.49751; // variational parameter
    const double step = 1.0; // only for metropolis
    const double dt = 0.1; // only for importance sampling
    // Mode 2 - varying alpha
    const double alpha_min = 0.3; // in mode 1 (varying alpha) minimum alpha
    const double alpha_max = 0.7; // in mode 2 (varying alpha) maximum alpha
    const int N_alpha = 5; // in mode 1 (varying alpha) number of different alphas between alpha_min and alpha_max
    const bool alpha_to_file = true; // set true to save data to file   

    // Select working mode : run "make" command. Then "./main x" where x is the number indicating working mode (see switch below)
    int selector = 0;
    if(argc>1){
        int a = stoi(argv[1]);
        assert(a>=0 && a<=6);
        selector = a;
    }

    // Information for the solvers
    const double D = 0.5; // only for importance sampling
    
    // Information for the hamiltonian
    const double a = 0.0043; // Set the radius of the particles. a=0 is the non-interacting case
    double omegaXY = 1.0; // Only the elliptical hamiltonian distinguish between omegaXY and omegaZ
    double omegaZ = sqrt(8); 

    // Information for the wavefunction
    const double beta = sqrt(8); // Only for asymmetrical wavefunction
    
    /*#######################################################################################################*/
    // Parameters for the various type of simulations
    // Mode 0 - normal

    // Mode 1 - numerical
    const double h = 1e-5; // Steplength for numerical derivatives and evaluations
    

    // Mode 3 - varying dt
    const double dt_min = 1e-3; // in mode 2 (varying dt) minimum dt
    const double dt_max = 10; // in mode 2 (varying dt) maximum dt
    const int N_dt = 20; // in mode 2 (varying dt) number of different dts between dt_min dt_max
    const bool dt_to_file = true; // set true to save data to file

    // Mode 4 - varying N
    vector<int> Ns {1, 10, 50, 100, 500}; // in mode 3 (varying N) different values of N
    const bool N_to_file = true; // set true to save data to file

    // Mode 5 - Gradient Descent
    double initial_alpha = 0.60;
    const double gamma = 1e-3; // Learning rate
    const unsigned int Nmax_gradient = 1000; // Max steps alowed
    const unsigned int Nsteps_gradient = (int) 1e5; // MC steps during the gradient descent
    const double tolerance = 1e-5; // Conditoin to stop the gradient descent

    // Mode 6 - One Body density
    const double r_max = 4.0;
    const int Nbins = 200;


    #if RUN_PARALLEL == 1
        int Nthreads = NTHREADS;
        int Ni = Nsteps_final/Nthreads;
    #else
        int Nthreads = 1;
        int Ni = Nsteps_final;
    #endif

    auto start = chrono::steady_clock::now(); // Store starting time to measure run time
    double x=0.0, y=0.0;

    omp_set_num_threads(Nthreads);
    #pragma omp parallel reduction (+:x,y)
    {
        System system((int) DIMENSION, (int) NPARTICLES, (bool) RUN_PARALLEL);

        // Hamiltonians
        #if USE_ELLIPTICAL == 1
            Elliptical elliptical(&system, omegaXY, omegaZ);
            system.setHamiltonian(&elliptical);
        #else
            Spherical spherical(&system, omegaXY);
            system.setHamiltonian(&spherical);
        #endif

        // Wavefunctions
        #if DIMENSION == 3 
            #if USE_ASYMMETRIC
                AsymmetricGaussian asymmgaussian(&system, alpha, beta, a);
                system.setWavefunction(&asymmgaussian);
            #else
                Gaussian gaussian(&system, alpha);
                system.setWavefunction(&gaussian);          
            #endif
        #else
            Gaussian gaussian(&system, alpha);
            system.setWavefunction(&gaussian);  
        #endif

        // Solvers
        #if USE_IMPORTANCE
            ImportanceSampling importance(&system, Ni, NstepsThermal, dt, D, (bool) TO_FILE);
            system.setSolver(&importance);
        #else
            Metropolis metropolis(&system, Ni, NstepsThermal, step, (bool) TO_FILE);
            system.setSolver(&metropolis);
        #endif

        // Others
        RandomGenerator randomgenerator;
        system.setRandomGenerator(&randomgenerator);

        Functions functions(&system, (bool) RUN_PARALLEL);
        
        if(omp_get_thread_num()==0) {functions.printPresentation();}
        
        vector<double> res(4,0.0);
        switch(selector){
            case 0: res=functions.solve_singleRun();  break; // Simple simulation
            case 1: functions.solve_singleRun(h); break; // Simple simulation with numerical derivative
            case 2: functions.solve_varying_alpha(alpha_min, alpha_max, N_alpha, alpha_to_file); break;
            case 3: functions.solve_varying_dt(dt_min, dt_max, N_dt, dt_to_file); break;
            case 4: functions.solve_varying_N(Ns, N_to_file); break;
            case 5: x+=functions.gradientDescent(initial_alpha, gamma, tolerance, Nmax_gradient, Nsteps_gradient); break;
            case 6: functions.solve_singleRun(r_max, Nbins); break;
        }
        x=res[0];
        y=res[2];

        
        

    }

    cout << "energy_averaged= " <<  x/Nthreads << "\t acceptance_averaged=" << y/Nthreads << endl;

    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << fixed << "Simulation termined. Total running time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl; // Print run time*/

    
    
}

