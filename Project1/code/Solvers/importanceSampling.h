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
        ~ImportanceSampling();
    
        // Getters
        double getdt();
        double getD();
    
        // Setters
        /// Set time steplength to integratenumerically the Langevin's equation
        void setdt(double dt);
        /// Set the diffusion coefficient in the Langevin's equation
        void setD(double D);
    
        // Other functions and attributes
        /// \see Solver::solve(bool AllAverages)
        vector<double> solve(bool allAverages); // uses LocalEnergyAnalytic()
        /// \see Solver::solve(double h)
        vector<double> solve(double h); // uses LocalEnergyNumeric()
    
    private:
        double dt;
        double D;
};