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
        double getStep();
        double setStep(double step);
        vector<double> solve(); // uses LocalEnergyAnalytic(), solves for a single alpha
        vector<double> solve(double h); // uses LocalEnergyNumeric(double h), solves for a single alpha
        vector<vector<double>> solve(double alpha_min, double alpha_max, int Nalphas); // uses LocalEnergyAnalytic(), solves for various alphas

    private:
        double step;
};