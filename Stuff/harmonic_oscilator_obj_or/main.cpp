#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <iomanip>

#include "system.h"
#include "gaussian.h"
#include "asymmetricGaussian.h"
#include "spherical.h"
#include "metropolis.h"
using namespace std;

int main(){

    auto start = chrono::steady_clock::now(); // Store starting time to measure run time
    srand(time(NULL));
    const double min = -1.0;
    const double max = 1.0;
    double alpha = 1;
    const double beta = 1.0;
    double omega = 1.0;
    const int dimension = 3;
    const int Nparticles = 100;

    const int Nsteps = (int) 1e4;
    const double step = 1.0;
    const double initialFraction = 0.1;

    System system(dimension, Nparticles);
    //Gaussian psi(&system, alpha); // Do we need this *here*?
    Spherical spherical(&system, omega);
    Metropolis metropolis(&system, Nsteps, step, initialFraction);

    //system.setWavefunction(&psi); // Do we need this *here*?
    system.setHamiltonian(&spherical);
    system.setSolver(&metropolis);



    // Different values of alpha
    double alpha_min = 0.2;
    double alpha_max = 0.8;
    int alpha_N = 10;

    
    // File to store results
    ofstream fout("results.csv");
    fout << "alpha,energy,STD" << endl;

    cout << endl << "Starting calculations... " << endl << endl;

    // Ground state energy estimation for each alpha, print data to results.csv
    for(int i=0; i<alpha_N; i++){
        alpha = alpha_min + i*(alpha_max - alpha_min)/alpha_N;
        system.setWavefunction(new Gaussian(&system, alpha));
        // remember that you can solve between analytical and numerical solver in metropolis.cpp
        // the numerical version evaluates numerically the second derivative that enters in the local energy definition.
        vector<double> estimations = system.solver->solve(); // Energy and standard deviation
        fout << alpha << "," << estimations[0] << "," << estimations[1] << endl;
        cout << setprecision(5) << "\t Cycle: " << i << "\t alpha: " << alpha << "\t energy: " << estimations[0] << "\t std: " << estimations[1] << endl;
    }

    fout.close();

    auto end = chrono::steady_clock::now(); // Store end time to measure run time
    auto diff = end - start; // Time difference
    cout << endl << "Simulation termined. Run time: " << chrono::duration <double, milli> (diff).count() << " ms" << endl << endl; // Print run time
    
    
    

}