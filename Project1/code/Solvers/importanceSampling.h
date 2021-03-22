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
        // Constructor and destructor
        ImportanceSampling(class System* system, int Nsteps, double initialFraction, double dt, double D, bool tofile);
        ~ImportanceSampling();
        
        // Other functions
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        vector<double> solve(double r_max, int N_bins); 
        vector<double> solve(double h); // uses LocalEnergyNumeric()
        
        void thermalize();
    
};