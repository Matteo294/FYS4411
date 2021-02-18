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
        virtual vector<double> computeEnergy() = 0;
        virtual vector<double> computeEnergy(double h) = 0; // override when numerical local energy
        virtual vector<double> computeExpectationPsidotPsiEl() = 0;
        virtual vector<double> computeExpectationPsidotPsi() = 0;
        class System* system;

    protected:
        int Nsteps;
        double InitialFraction;

};