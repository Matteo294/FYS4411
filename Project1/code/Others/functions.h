#pragma once
#include "../System/system.h"
/*#include "../System/system.h"
#include "../Wavefunctions/gaussian.h"
#include "../Wavefunctions/asymmetricGaussian.h"
#include "../Hamiltonians/spherical.h"
#include "../Hamiltonians/elliptical.h"
#include "../Solvers/metropolis.h"
#include "../Solvers/importanceSampling.h"
#include "../Others/random_generator.h"*/
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <cmath>
#include <experimental/filesystem>

class Functions{

    public:
        Functions(class System* system, bool parallel);
        ~Functions();

        class System* system;

        vector<double> solve_singleRun();
        vector<double> solve_singleRun(double h);
        vector<double> solve_singleRun(double r_max, int Nbins);


        // functions
        /** This function calls recursively the solver for Nalphas different alpha values between alpha_min and alpha_max. 
         * Remember that the parameter alpha is the variational parameter of the gaussian wavefunction
        **/
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool alphatoFile); 
        
        /** This function calls recursively the solver (ImportanceSampling only!) for Ndt different values of dt between dt_min and dt_max.
         * Remember that the parameter dt is the one used to integrate the Langevin's equation numerically
         * \attention This function works only with the importance sampling solver
         * \see ImportanceSampling::ImportanceSampling()
        **/
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt, bool dttoFile);

        /** This function calls recursively the solver for the given values of number of particles (pass the single values via vector)
        **/
        vector<vector<double>> solve_varying_N(vector<int> N, bool NtoFile);       

        /// Finds the best alpha parameter using gradient descent 
        double gradientDescent(double initialAlpha, double gamma, double tolerance, int Nmax, int Nsteps);

        void printPresentation();

        /// Useful function to print in a standard way the results given by the solver
        void printResultsSolver(vector<double> res);
        
    private:
        ofstream alphaFile, dtFile, NFile;
        bool parallel;

};