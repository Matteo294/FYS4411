#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

class ImportanceSampling : public Solver{
    public:
        ImportanceSampling(class System* system, int Nsteps, double initialFraction, double dt, double D);
    
    // Getters
        double getdt();
        double getD();
    
    // Setters
        void setdt(double dt);
        void setD(double D);
    
    // Other functions and attributes
        vector<double> computeEnergy(); // uses LocalEnergyAnalytic()
        vector<double> computeEnergy(double h); // uses LocalEnergyNumeric()
        vector<double> computeExpectationPsidotPsi(); // computes average psi'/psi -> can be made more efficient by using logarithm
        vector<double> computeExpectationPsidotPsiEl(); // computes average psi'/psi*Eloc -> can be made more efficient by using logarithm
    
    private:
        double dt;
        double D;
};