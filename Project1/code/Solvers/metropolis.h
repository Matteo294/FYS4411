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
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        vector<double> solve(double h); // uses LocalEnergyNumeric(double h)

};