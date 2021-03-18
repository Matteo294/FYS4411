#include "config.h"
using namespace std;

// Parallel or serial version of the code (check compatibility with the chosen working mode)
bool parallelize;

// Select pieces of the system: wavefunction (symmetric-asymetric), hamiltonian (spherical-elliptical), solver(bruteforce-metropolis)
bool activate_asymmetric;
bool activate_elliptical;
bool activate_importance;

// Information for the system
const int dimension = 3;
const int Nparticles = 5;

// Information for the solvers
const int Nsteps_final = (int) pow(2,21); // MC steps for the final simulation
const int NstepsThermal = (int) 1e5; // Fraction of septs to wait for the system thermalization
const double step = 1.0; // only for metropolis
const double D = 0.5; // only for importance sampling
const double dt = 0.01; // only for importance sampling

// Information for the hamiltonian
const double a = 0.043; // Set the radius of the particles. a=0 is the non-interacting case
double omegaXY = 1.0; // Only the elliptical hamiltonian distinguish between omegaXY and omegaZ
double omegaZ = 2.82843; 

// Information for the wavefunction
double alpha = 0.5; // variational parameter
const double beta = 2.82843; // Only for asymmetrical wavefunction

// Others
bool tofile = true; // Print on external file for resampling analysis (numerical methods)

// Parameters for the various type of simulations

// Mode 0 - normal

// Mode 1 - numerical
const double h = 1e-5; // Steplength for numerical derivatives and evaluations

// Mode 2 - varying alpha
const double alpha_min = 0.3; // in mode 1 (varying alpha) minimum alpha
const double alpha_max = 0.7; // in mode 2 (varying alpha) maximum alpha
const int N_alpha = 20; // in mode 1 (varying alpha) number of different alphas between alpha_min and alpha_max
const bool alpha_to_file = false; // set true to save data to file

// Mode 3 - varying dt
const double dt_min = 1e-3; // in mode 2 (varying dt) minimum dt
const double dt_max = 10; // in mode 2 (varying dt) maximum dt
const int N_dt = 50; // in mode 2 (varying dt) number of different dts between dt_min dt_max
const bool dt_to_file = false; // set true to save data to file

// Mode 4 - varying N
vector<int> Ns {5, 10, 15}; // in mode 3 (varying N) different values of N
const bool N_to_file = false; // set true to save data to file

// Mode 5 - Gradient Descent
double best_alpha = 0.0;
double initial_alpha = 0.49;
//const double gamma = 2e-3; // Learning rate
const unsigned int Nmax_gradient = 100; // Max steps alowed
const unsigned int Nsteps_gradient = (int) 1e5; // MC steps during the gradient descent
const double tolerance = 1e-8; // Conditoin to stop the gradient descent

// Mode 6 - One Body density
const double r_max = 4.0;
int Nbins = 200;