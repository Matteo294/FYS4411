#pragma once 
#include "wavefunction.h"

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha); // Constructor
        double evaluate(class System* s);
        double evaluateSecondDerivative(class System* s);
        double numericalSecondDerivative();
    private:
        double alpha; // Variational parameter
};