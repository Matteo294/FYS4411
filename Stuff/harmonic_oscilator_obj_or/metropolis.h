#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

class Metropolis : public Solver{
    public:
        Metropolis(class System* system, int Nsteps, double step, double initialFraction);
        int getNsteps();
        double getstep();
        double getinitialFraction();
        vector<double> solve();
        double best_alpha;

};