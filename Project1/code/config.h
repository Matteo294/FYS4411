#pragma once 
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
#include <string>
#include <cmath>

// Parallel or serial version of the code (check compatibility with the chosen working mode)
extern bool parallelize;

// Select pieces of the system: wavefunction (symmetric-asymetric), hamiltonian (spherical-elliptical), solver(bruteforce-metropolis)
extern bool activate_asymmetric;
extern bool activate_elliptical;
extern bool activate_importance;

// Information for the system
extern const int dimension;
extern const int Nparticles;

// Information for the solvers
extern const int Nsteps_final; // MC steps for the final simulation
extern const int NstepsThermal; // Fraction of septs to wait for the system thermalization
extern const double step; // only for metropolis
extern const double D; // only for importance sampling
extern const double dt; // only for importance sampling

// Information for the hamiltonian
extern const double a; // Set the radius of the particles. a=0 is the non-interacting case
extern double omegaXY; // Only the elliptical hamiltonian distinguish between omegaXY and omegaZ
extern double omegaZ; 

// Information for the wavefunction
extern double alpha; // variational parameter
extern const double beta; // Only for asymmetrical wavefunction

// Others
extern bool tofile; // Print on external file for resampling analysis (numerical methods)

// Parameters for the various type of simulations

// Mode 0 - normal

// Mode 1 - numerical
extern const double h; // Steplength for numerical derivatives and evaluations

// Mode 2 - varying alpha
extern const double alpha_min; // in mode 1 (varying alpha) minimum alpha
extern const double alpha_max; // in mode 2 (varying alpha) maximum alpha
extern const int N_alpha; // in mode 1 (varying alpha) number of different alphas between alpha_min and alpha_max
extern const bool alpha_to_file; // set true to save data to file

// Mode 3 - varying dt
extern const double dt_min; // in mode 2 (varying dt) minimum dt
extern const double dt_max; // in mode 2 (varying dt) maximum dt
extern const int N_dt; // in mode 2 (varying dt) number of different dts between dt_min dt_max
extern const bool dt_to_file; // set true to save data to file

// Mode 4 - varying N
extern vector<int> Ns; // in mode 3 (varying N) different values of N
extern const bool N_to_file; // set true to save data to file

// Mode 5 - Gradient Descent
extern double best_alpha;
extern double initial_alpha;
//extern const double gamma; // Learning rate
extern const unsigned int Nmax_gradient; // Max steps alowed
extern const unsigned int Nsteps_gradient; // MC steps during the gradient descent
extern const double tolerance; // Conditoin to stop the gradient descent

// Mode 6 - One Body density
extern const double r_max;
extern int Nbins;