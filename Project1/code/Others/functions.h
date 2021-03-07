#pragma once
#include "../System/system.h"
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>

class Functions{

    public:
        Functions(class System* system);
        ~Functions();

        class System* system;

        // functions
        /** This function calls recursively the solver (set it in the system class via the proper setter) for Nalphas different alpha values
         * betwee alpha_min and alpha_max. Remember that the parameter alpha is the variational parameter of the gaussian wavefunction
        **/
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool toFile); 
        
        /** This function calls recursively the solver (ImportanceSampling only!!!) for Ndt different values of dt between dt_min and dt_max.
         * Remember that the parameter dt is the one used to integrate the Langevin's equation numerically
         * \attention This function works only with the importance sampling solver
         * \see ImportanceSampling::ImportanceSampling()
        **/
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt); 

        /** This function calls recursively the solver (set it in the system class via the proper setter) for the values of number of particles
         * given via the vector N.
         * \note The parameters (for example those of the solver and wavefunction) should be set in advance
        **/
        vector<vector<double>> solve_varying_N(vector<int> N);

        /// Finds the best alpha parameter using gradient descent 
        double gradientDescent(double initialAlpha, double gamma, double tolerance, int Nmax, int Nsteps);

        /** Prints the value to the file followed by the string s. Set c as the separation character (e.g. "," for csv) or newline character
         * \note That the file should already be open
        **/
        void printToFile(double val, ofstream f, string s);
    private:
        ofstream alphaFile, dtFile, Nfile;

};