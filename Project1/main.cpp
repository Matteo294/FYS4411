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
#include "Solvers/metropolis.h"
#include "Solvers/importanceSampling.h"
#include "Others/random_generator.h"
#include "Others/functions.h"
using namespace std;

int main(){

    double omega = 1.0;
    const int dimension = 3;
    const int Nparticles = 10;

    const int Nsteps = (int) 1e6;
    const double initialFraction = 0.1;
    const double step = 1.0; // only for metropolis
    const double D = 0.5; // only for importance sampling
    const double dt = 0.001; // only for importance sampling
    const double alpha = 0.5;
    

    System system(dimension, Nparticles);
    Spherical spherical(&system, omega);
    Gaussian gaussian(&system, alpha);
    AsymmetricGaussian wf(&system, alpha, (double) 5.0);
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
    
    double initialAlpha = 0.4;
    double gamma = 1e-2;
    double tolerance = 1e-8;
    int Nmax = 50;
    int Nsteps_gradient = (int) 1e4;
    functions.gradientDescent(initialAlpha, gamma, tolerance, Nmax, Nsteps_gradient);

    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << "Simulation termined. Simulation time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl << endl; // Print run time


}