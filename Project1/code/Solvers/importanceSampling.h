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
        ImportanceSampling(class System* system, int Nsteps, double initialFraction, double dt, double D, bool tofile);
        ~ImportanceSampling();
        
        vector<double> solve(bool allAverages);
        vector<double> solve(double r_max, int N_bins); 
        vector<double> solve(double h); 
        
        void thermalize();
    
};