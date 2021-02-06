#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
<<<<<<< HEAD
=======
#include <iomanip>
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
using namespace std;

class ImportanceSampling : public Solver{
    public:
        ImportanceSampling(class System* system, int Nsteps, double initialFraction, double dt, double D);
    
    // Getters
        double getdt();
        double getD();
    
    // Setters
        void setdt(double dt);
        void setD(double D);
    
    // Other functions and attributes
        vector<double> solve(); // uses LocalEnergyAnalytic()
        vector<double> solve(double h); // uses LocalEnergyNumeric()
<<<<<<< HEAD
=======
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas);
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt);
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
    
    private:
        double dt;
        double D;
<<<<<<< HEAD


=======
>>>>>>> 6c76b22c10e7e8f3688d66d11aab33e7d1666a61
};