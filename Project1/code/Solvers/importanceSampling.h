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
    
    // Other functions and attributes
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        vector<double> solve(double h); // uses LocalEnergyNumeric()
    
};