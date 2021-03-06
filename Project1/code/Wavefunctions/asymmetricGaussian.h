#pragma once
#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta, double a); // Constructor
        double evaluateAll();
        double evaluateSing(int part_idx);
        double psibar_psi();
        double numericalSecondDerivative(int part_idx, int direction, double h);
        vector<double> DriftForce(int part_idx);
};