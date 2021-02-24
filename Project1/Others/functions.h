#pragma once
#include "../System/system.h"
#include <vector>
#include <iomanip>

class Functions{

    public:
        // creator
        Functions(class System* system);
        class System* system;

        // functions
        // solves for a bunch of alpha values, uses solve((bool) 0) --> computes only energy and std for energy
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas); 
        
        // solves for a bunch of dt values, mandatory having chosen ImportanceSampling solver
        //vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt); 

        // finds the best alpha parameter using gradient descent 
        double gradientDescent(double inizialAlpha, double gamma, double tolerance, int Nmax, int Nsteps);


};