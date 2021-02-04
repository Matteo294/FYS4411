#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
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
    
    private:
        double dt;
        double D;


};