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
    double alpha = 0.60; // variational parameter
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

    auto start = chrono::steady_clock::now(); // Store starting time to measure run time
    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference


    for(int n=2; n<100; n++){
        cout << "N particles: " << n << "\t";
        system.addParticle(1.0, {0.0, 0.0, 0.0});

        // Gradient descent
        start = chrono::steady_clock::now(); // Store starting time to measure run time
        alpha = functions.gradientDescent(alpha, gamma, tolerance, Nmax_gradient, Nsteps);
        stop = chrono::steady_clock::now(); // Store starting time to measure run time
        diff = stop - start; // Time difference
        cout << "Best alpha found. Time taken: " << chrono::duration <double, milli> (diff).count()/1000 << " s \t Best alpha: " << alpha << endl; // Print run time

        // Final simulation
        system.getSolver()->setNsteps(Nsteps_final);
        start = chrono::steady_clock::now();
        vector<double> res = system.getSolver()->solve(false);
        stop = chrono::steady_clock::now();
        cout << fixed << setprecision(5) << res[0] << "\t" << res[1] << "\t" << res[2] << endl;
        diff = stop - start; // Time difference
        
        //functions.solve_varying_alpha((double) 0.3, (double) 0.7, (int) 4);
        
        cout << endl << "Simulation termined. Simulation time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl; // Print run time

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

}