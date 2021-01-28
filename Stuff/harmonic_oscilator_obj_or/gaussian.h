#pragma once 
#include "wavefunction.h"
#include <cmath>

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha); // Constructor
        double evaluate();
        double evaluateSecondDerivative();
        double numericalSecondDerivative();
        double alpha; // Variational parameter
};