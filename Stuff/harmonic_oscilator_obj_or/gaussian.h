#pragma once 
#include "wavefunction.h"
#include <cmath>

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha); // Constructor
        virtual double evaluate();
        virtual double evaluateSecondDerivative();
        virtual double numericalSecondDerivative();
    private:
        double alpha; // Variational parameter
};