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
        vector<double> computeEnergy(); // uses LocalEnergyAnalytic()
        vector<double> computeEnergy(double h); // uses LocalEnergyNumeric(double h)

    private:
        double step;
};