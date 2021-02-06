#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
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
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas);
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt);
    
    private:
        double dt;
        double D;
};