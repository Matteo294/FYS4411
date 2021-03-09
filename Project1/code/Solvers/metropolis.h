#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

class Metropolis : public Solver{
    public:
        Metropolis(class System* system, int Nsteps, double initialFraction, double step);
        ~Metropolis();

        // Getters and setters
        double getStep();
        double setStep(double step);

        /// \see Solver::solve(bool allAverages)
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        /// \see Solver::solve(double r_max, int N_bins)
        vector<double> solve(double r_max, int N_bins);
        /// \see Solver::solve(double h)
        vector<double> solve(double h); // uses LocalEnergyNumeric(double h)


};