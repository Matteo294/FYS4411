#pragma once
#include "../System/system.h"
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

        /** Runs a single Monte Carlo simulation using the analytic expression for the local energy.
         * \see Solver::solve() **/
        vector<double> solve_singleRun();
        /** Runs a single Monte Carlo simulation using numerical evaluation of the second derivative
         * appearing in the local energy.
         * \see Solver::solve(double h) **/
        vector<double> solve_singleRun(double h);
        /** Runs a single Monte Carlo simulation using the analytic expression for the local energy and producing the histogram for the
         * radial one-body density. 
         * \see Solver::solve(double r_max, int Nbins) **/
        vector<double> solve_singleRun(double r_max, int Nbins);


        // functions
        /** This function calls recursively the solver for Nalphas+1 different alpha values between alpha_min and alpha_max. 
         * Remember that the parameter alpha is the variational parameter of the wavefunction.
         * Within the function the user can set the path where to generate the two .csv files containing respectively the counts for each bin
         * and the alpha values.
         * \attention The alpha values are printed to a .csv only if bool alphatofile=True.
         * \note The .csv file with the counts for each bin is always generated and filled.
        **/
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool alphatoFile); 
        
        /** This function calls recursively the solver (ImportanceSampling only!) for Ndt+1 different values of dt between dt_min and dt_max.
         * Remember that the parameter dt is the one used to integrate the Langevin's equation numerically.
         * Within the function the user can set the path where to generate the .csv file containing the values of dt and also the
         * .dat files with the local energy values saved at every step for each dt value.
         * \attention The dt values are printed to a .csv only if bool dttofile=True.
         * \note The .dat files with the local energy value at every step are generated only if the toFile flag in Solver is active.
         * \see Solver
         * \attention This function works only with the importance sampling solver.
         * \see ImportanceSampling::ImportanceSampling()
        **/
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt, bool dttoFile);

        /** This function progressively adds particles to the system for the values specified in the argument vector and
         * calls recursively the solver.
         * Within the function the user can set the path where to generate the .csv file containing the values of N and also the
         * .dat files with the local energy values saved at every step for each N value.
         * \attention The N values are printed to a .csv only if bool Ntofile=True.
         * \note The .dat files with the local energy value at every step are generated only if the toFile flag in Solver is active.
         * \see Solver
         * \attention This function works only if the N values specified in the vector are increasing and if the minimum among them
         * (the first in the vector) is smaller or equal to the number of particle that we set in the system at the moment of its creation.
         * \see System::system()
        **/
        vector<vector<double>> solve_varying_N(vector<int> N, bool NtoFile);       

        /** Finds the best alpha parameter using gradient descent approach starting from an initialAlpha. 
         * gamma specifies the learning rate
         * tolerance specifies the condition on the derivative to stop the execution
         * Nmax is the maximum number of iterations allowed before interrupting the execution
         * Nsteps is the number of steps used in each Monte Carlo simulation.
         * \see Solver::solve(true) **/
        double gradientDescent(double initialAlpha, double gamma, double tolerance, int Nmax, int Nsteps);

        /// Print names of the authors of the code at the beginning of every execution
        void printPresentation();

        /** Useful function to print in a standard way the results given by the solver
         * \see Solver::solve() **/        
        void printResultsSolver(vector<double> res);
        
        /// Print the configuration and the parameters of the system before launching the simulation
        void printConfiguration(int selector, bool asymmetric, bool elliptical, bool importance);

    private:
        ofstream alphaFile, dtFile, NFile;
        bool parallel;

};