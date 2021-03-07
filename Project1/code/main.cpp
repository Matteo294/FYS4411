#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <iomanip>

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

    // Select working mode
    int selector = 0;
    if(argc>1){
        int a = stoi(argv[1]);
        assert(a>0 && a<2);
        selector = stoi(argv[1]);
    }

    // Information for the system
    const int dimension = 3;
    const int Nparticles = 1;

    // Information for the solvers
    const int Nsteps_final = (int) 1e6; // MC steps for the final simulation
    const double step = 1.0; // only for metropolis
    const double D = 0.5; // only for importance sampling
    const double dt = 0.01; // only for importance sampling
    
    // Information for the hamiltonian
    const double a = 0.0; // Set the radius of the particles. a=0 is the non-interacting case
    double omegaXY = 1.0; // Only the elliptical hamiltonian distinguish between omegaXY and omegaZ
    double omegaZ = 1.0; 

    // Information for the wavefunction
    double alpha = 0.50; // variational parameter
    const double beta = 1.0; // Only for asymmetrical wavefunction

    // Gradient descent
    const double gamma = 1e-2; // Learning rate
    const unsigned int Nmax_gradient = 100; // Max steps alowed
    const unsigned int Nsteps = (int) 1e4; // MC steps during the gradient descent
    const double tolerance = 1e-8; // Conditoin to stop the gradient descent
    
    // Others
    const double h = 1e-5; // Steplength for numerical derivatives and evaluations
    const float initialFraction = 0.1; // Fraction of septs to wait for the system thermalization
    bool dt_analysis = false; // Plotting flags: turn True to save data to make the plots

    // Parameters for the various type of simulations
    // Mode 1 - varying alpha
    const double alpha_min = 0.2; // in mode 1 (varying alpha) minimum alpha
    const double alpha_max = 0.8; // in mode 2 (varying alpha) maximum alpha
    const int N_alpha = 100; // in mode 1 (varying alpha) number of different alphas between alpha_min and alpha_max
    const bool alpha_to_file = true; // set true to save data to file
    // Mode 2 - varying dt
    const double dt_min = 1e-8; // in mode 2 (varying dt) minimum dt
    const double dt_max = 1e-2; // in mode 2 (varying dt) maximum dt
    const int N_dt = 10; // in mode 2 (varying dt) number of different dts between dt_min dt_max
    const bool dt_to_file = false; // set true to save data to file
    // Mode 3 - varying N
    vector<int> Ns {2, 5, 10}; // in mode 3 (varying N) different values of N
    const bool N_to_file = false; // set true to save data to file

    System system(dimension, Nparticles);

    // Hamiltonians
    Spherical spherical(&system, omegaXY);
    Elliptical elliptical(&system, omegaXY, omegaZ);

    // Wavefunctions
    Gaussian gaussian(&system, alpha);
    AsymmetricGaussian asymmgaussian(&system, alpha, beta, a);

    // Solvers
    Metropolis metropolis(&system, Nsteps, initialFraction, step);
    ImportanceSampling importance(&system, Nsteps, initialFraction, dt, D);
    
    // Others
    RandomGenerator randomgenerator;
    Functions functions(&system);

    // Choose options
    system.setHamiltonian(&elliptical);
    system.setWavefunction(&asymmgaussian);
    system.setSolver(&metropolis);
    system.setRandomGenerator(&randomgenerator);

    functions.printPresentation();

    // !!!!!!!! This should be more precise and should be put inside Function (time should be returned as another value)
    auto start = chrono::steady_clock::now(); // Store starting time to measure run time

    switch(selector){
        case 0: functions.printResultsSolver(system.getSolver()->solve(false)); break; // Simple simulation
        case 1: functions.solve_varying_alpha(alpha_min, alpha_max, N_alpha, alpha_to_file); break;
        case 2: functions.solve_varying_dt(dt_min, dt_max, N_dt, dt_to_file); break;
        case 3: functions.solve_varying_N(Ns, N_to_file); break;
    }
    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << "Simulation termined. Total running time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl; // Print run time*/

    // Sasha's part
    /*ofstream myfile;
    myfile.open ("energy.dat");
    //myfile << "energy" << endl;
    int i; 
    for(i=0;i<100;i++){
        vector<double> res = system.getSolver()->solve(h);
        myfile << fixed << setprecision(24) << res[0] << endl;
        cout << fixed << setprecision(24) << res[0] << "\t" << res[1] << "\t" << res[2] << endl;
    }*/

}