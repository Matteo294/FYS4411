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
<<<<<<< HEAD
    double alpha = 1.0;
    const double beta = 1.0;
=======
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
    double omega = 1.0;
    const int dimension = 3;
    const int Nparticles = 10;

<<<<<<< HEAD
    const int Nsteps = (int) 1e5;
=======
    const int Nsteps = (int) 1e4;
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
    const double step = 1.0;
    const double initialFraction = 0.1;
    const double D = 0.5;
    const double dt = 0.01;

    System system(dimension, Nparticles);
    Spherical spherical(&system, omega);
<<<<<<< HEAD
=======
    Gaussian gaussian(&system, (double) 0.5);
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
    //Metropolis metropolis(&system, Nsteps, initialFraction, step);
    ImportanceSampling importance(&system, Nsteps, initialFraction, dt, D);
    RandomGenerator randomgenerator;

    system.setHamiltonian(&spherical);
<<<<<<< HEAD
=======
    system.setWavefunction(&gaussian);
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
    //system.setSolver(&metropolis);
    system.setSolver(&importance);
    system.setRandomGenerator(&randomgenerator);


    // Different values of alpha
    double alpha_min = 0.2;
    double alpha_max = 0.8;
    int alpha_N = 10;

<<<<<<< HEAD
    
    // File to store results
    ofstream fout("results.csv");
    fout << "alpha,energy,STD" << endl;

    cout << endl << "Starting calculations... " << endl << endl;

    // Ground state energy estimation for each alpha, print data to results.csv
    for(int i=0; i<alpha_N; i++){
        alpha = alpha_min + i*(alpha_max - alpha_min)/alpha_N;
        system.setWavefunction(new Gaussian(&system, alpha));
        vector<double> estimations = system.getSolver()->solve(); // Energy and standard deviation
        fout << alpha << "," << estimations[0] << "," << estimations[1] << endl;
        cout << setprecision(5) << "\t Cycle: " << i << "\t alpha: " << alpha << "\t energy: " << estimations[0] << "\t std: " << estimations[1] << endl;
    }

    fout.close();

    auto end = chrono::steady_clock::now(); // Store end time to measure run time
    auto diff = end - start; // Time difference
    cout << endl << "Simulation termined. Run time: " << chrono::duration <double, milli> (diff).count() << " ms" << endl << endl; // Print run time
    
    
=======
    // different values of dt
    double dt_min = 0.001;
    double dt_max = 0.1;
    int Ndt = 10;

    vector<vector<double>> results;
    results = system.getSolver()->solve_varying_dt(dt_min, dt_max, Ndt);   
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
    

}