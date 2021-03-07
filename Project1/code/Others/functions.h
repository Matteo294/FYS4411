#pragma once
#include "../System/system.h"
#include <vector>
#include <iomanip>
#include <cmath>

class Functions{

    public:
        Functions(class System* system);
        ~Functions();

        class System* system;

        // functions
        /** The function calls recursively the solver (set it in the system class via the proper setter) for Nalphas different alpha values
         * betwee alpha_min and alpha_max. Remember that the parameter alpha is the variational parameter of the gaussian wavefunction
        **/
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas); 
        
        /** The function calls recursively the solver (ImportanceSampling only!!!) for Ndt different values of dt between dt_min and dt_max.
         * Remember that the parameter dt is the one used to integrate the Langevin's equation numerically
         * \attention This function works only with the importance sampling solver
         * \see ImportanceSampling::ImportanceSampling()
        **/
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt); 

        // finds the best alpha parameter using gradient descent 
        double gradientDescent(double inizialAlpha, double gamma, double tolerance, int Nmax, int Nsteps);


};