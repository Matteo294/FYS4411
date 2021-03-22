#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

class Metropolis : public Solver{
    public:
        // Constructor and destructor
        Metropolis(class System* system, int Nsteps, int NstepsThermal, double step, bool tofile);
        ~Metropolis();

        // Other functions
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        vector<double> solve(double r_max, int N_bins);
        vector<double> solve(double h); // uses LocalEnergyNumeric(double h)

        void thermalize();


};