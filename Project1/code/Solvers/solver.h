#pragma once
#include "../System/system.h"
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
        // performs a MC simulation to evaluate the energy and the std on the energy.
        // if allAverages = 1 --> evaluates the derivative of the mean value of the energy with respect to alpha
        virtual vector<double> solve(bool allAverages) = 0;
        virtual vector<double> solve(double h) = 0; // override when numerical local energy
        class System* system;

    protected:
        int Nsteps;
        double InitialFraction;

};