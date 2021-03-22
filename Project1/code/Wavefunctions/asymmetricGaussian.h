#pragma once
#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta, double a);
        ~AsymmetricGaussian();
        
        double evaluateAll();
        double evaluateSing(int part_idx);
        double numericalSecondDerivative(int part_idx, int direction, double h);
        double psibar_psi();
        vector<double> DriftForce(int part_idx);
};