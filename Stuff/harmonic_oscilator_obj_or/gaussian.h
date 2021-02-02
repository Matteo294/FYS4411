#pragma once 
#include "wavefunction.h"
#include <cmath>

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha); // Constructor
        double evaluateAll();
        double evaluateSing(int part_idx);
        double evaluateSecondDerivative();
        double numericalSecondDerivative(int part_idx, int direction, double h);
        virtual
        vector<double> DriftForce(int part_idx);
        double alpha; // Variational parameter
};