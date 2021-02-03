#pragma once 
#include "wavefunction.h"
#include <cmath>
#include <vector>

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha); // Constructor
        double evaluateAll();
        double evaluateSing(int part_idx);
        double evaluateSecondDerivative();
        std::vector<double> evaluateGradient();
        double numericalGradient(int part_idx, int direction, double h);
        double numericalSecondDerivative(int part_idx, int direction, double h);
        double alpha; // Variational parameter
};