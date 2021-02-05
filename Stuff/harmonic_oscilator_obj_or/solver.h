#pragma once
#include "system.h"
#include <ctime>
#include <vector>
using namespace std;

class Solver{
    public:
        Solver(class System* system, int Nsteps, double initialFraction);

    // Getters
        int getNsteps();
        double getInitialFraction();

    // Setters
        void setNsteps(int Nsteps);
        void setInitialFraction(double initialFraction);

    // Other functions & attributes
        virtual vector<double> solve() = 0; // solves for a single alpha
        virtual vector<double> solve(double h) = 0; // override when numerical local energy
        virtual vector<vector<double>> solve(double alpha_min, double alpha_max, int Nalphas) = 0; // solves for a bunch of alpha values
        class System* system;

    protected:
        int Nsteps;
        double InitialFraction;

};