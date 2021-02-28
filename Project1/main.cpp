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

int main(){

    double omegaXY = 1.0;
    double omegaZ = 1.3;
    const int dimension = 3;
    const int Nparticles = 10;

    const int Nsteps = (int) 1e4;
    const int Nsteps_final = (int) 1e7;
    const double initialFraction = 0.1;
    const double step = 1.0; // only for metropolis
    const double D = 0.5; // only for importance sampling
    const double dt = 0.01; // only for importance sampling
    const double alpha = 0.5;
    const double beta = 0.7;
    const double a = 0.0;
    
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
    system.setHamiltonian(&spherical);
    system.setWavefunction(&gaussian);
    system.setSolver(&metropolis);
    //system.setSolver(&importance);
    system.setRandomGenerator(&randomgenerator);

    auto start = chrono::steady_clock::now(); // Store starting time to measure run time
    
    functions.gradientDescent((double) 0.4, (double) 0.01, (double) 1e-8,  (int) 50, (int) 1e4);

    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << "Simulation termined. Simulation time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl << endl; // Print run time


}