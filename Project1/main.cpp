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
using namespace std;

int main(){

    auto start = chrono::steady_clock::now(); // Store starting time to measure run time
    double omega = 1.0;
    const int dimension = 3;
    const int Nparticles = 10;

    const int Nsteps = (int) 1e4;
    const double step = 1.0;
    const double initialFraction = 0.1;
    const double D = 0.5;
    const double dt = 0.01;

    System system(dimension, Nparticles);
    Spherical spherical(&system, omega);
    Gaussian gaussian(&system, (double) 0.5);
    //Metropolis metropolis(&system, Nsteps, initialFraction, step);
    ImportanceSampling importance(&system, Nsteps, initialFraction, dt, D);
    RandomGenerator randomgenerator;

    system.setHamiltonian(&spherical);
    system.setWavefunction(&gaussian);
    //system.setSolver(&metropolis);
    system.setSolver(&importance);
    system.setRandomGenerator(&randomgenerator);


    // Different values of alpha
    double alpha_min = 0.2;
    double alpha_max = 0.8;
    int alpha_N = 10;

    // different values of dt
    double dt_min = 0.001;
    double dt_max = 0.1;
    int Ndt = 10;

    vector<vector<double>> results;
    results = system.getSolver()->solve_varying_dt(dt_min, dt_max, Ndt);   
    

}