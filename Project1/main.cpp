#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <iomanip>

#include "System/system.h"
#include "Wavefunctions/gaussian.h"
//#include "Wavefunctions/asymmetricGaussian.h"
#include "Hamiltonians/spherical.h"
#include "Solvers/metropolis.h"
#include "Solvers/importanceSampling.h"
#include "Others/random_generator.h"
#include "Others/functions.h"
using namespace std;

int main(){

    double omega = 1.0;
    const int dimension = 1;
    const int Nparticles = 1;

    const int Nsteps = (int) 1e3;
    const double step = 1.0;
    const double initialFraction = 0.1;
    const double D = 0.5;
    const double dt = 0.001;
    const double alpha = 0.5;
    

    System system(dimension, Nparticles);
    Spherical spherical(&system, omega);
    Gaussian gaussian(&system, alpha);
    //Metropolis metropolis(&system, Nsteps, initialFraction, step);
    ImportanceSampling importance(&system, Nsteps, initialFraction, dt, D);
    RandomGenerator randomgenerator;
    Functions functions(&system);

    system.setHamiltonian(&spherical);
    system.setWavefunction(&gaussian);
    //system.setSolver(&metropolis);
    system.setSolver(&importance);
    system.setRandomGenerator(&randomgenerator);

    auto start = chrono::steady_clock::now(); // Store starting time to measure run time
    
    system.getSolver()->solve(); 
    
    // Different values of alpha
    double alpha_min = 0.2;
    double alpha_max = 0.8;
    int Nalphas = 10;
    
    // different values of dt
    /*double dt_min = 0.001;
    double dt_max = 0.1;
    int Ndt = 10; */


    

    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    functions.solve_varying_alpha(alpha_min, alpha_max, Nalphas);
    auto diff = stop - start; // Time difference
    cout << endl << "Simulation termined. Simulation time: " << chrono::duration <double, milli> (diff).count() << " ms" << endl << endl; // Print run time


}