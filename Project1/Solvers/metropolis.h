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
        vector<double> solve(); // uses LocalEnergyAnalytic()
        vector<double> solve(double h); // uses LocalEnergyNumeric(double h)
        vector<vector<double>> solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas); // uses LocalEnergyAnalytic(), solves for various alphas
        vector<vector<double>> solve_varying_dt(double dt_min, double dt_max, int Ndt);

    private:
        double step;
};