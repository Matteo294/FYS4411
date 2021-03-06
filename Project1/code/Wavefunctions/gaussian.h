#pragma once 
#include "wavefunction.h"
#include <cmath>

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha);
        ~Gaussian();
        
        double evaluateAll();
        double evaluateSing(int part_idx);
        double numericalSecondDerivative(int part_idx, int direction, double h);
        double psibar_psi();
        vector<double> DriftForce(int part_idx);
};