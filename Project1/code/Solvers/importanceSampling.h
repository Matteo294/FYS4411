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
        ~ImportanceSampling();
        
        // Other functions and attributes
        /// \see Solver::solve(bool AllAverages)
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        /// \see Solver::solve(double r_max, int N_bins)
        vector<double> solve(double r_max, int N_bins); 
        /// \see Solver::solve(double h)
        vector<double> solve(double h, bool tofile); // uses LocalEnergyNumeric()
    
};