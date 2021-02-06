#pragma once
#include "../System/system.h"
#include <ctime>
#include <vector>
using namespace std;

class Solver{
    public:
        Solver(class System* system, int Nsteps, double initialFraction);

    // Getters
        int getNsteps();
        double getInitialFraction();

    // Setters
        void setNsteps(int Nsteps);
        void setInitialFraction(double initialFraction);

    // Other functions & attributes
        virtual vector<double> solve() = 0;
        virtual vector<double> solve(double h) = 0; // override when numerical local energy
<<<<<<< HEAD
=======
        virtual vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas) = 0; // solves for a bunch of alpha values
        virtual vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt) = 0; 
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
        class System* system;

    protected:
        int Nsteps;
        double InitialFraction;

};