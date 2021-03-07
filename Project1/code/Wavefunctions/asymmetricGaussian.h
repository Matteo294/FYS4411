#pragma once
#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta, double a);
        ~AsymmetricGaussian();
        
        /// \see Wavefunction::evaluateAll()
        double evaluateAll();
        /// \see Wavefunctino::evaluateSing(int part_idx)
        double evaluateSing(int part_idx);
        /// \see Wavefunction::analyticalAlphaDerivative()
        double analyticalAlphaDerivative();
        /// \see Wavefunction::numericalSecondDerivative(int part_idx, int direction, double h)
        double numericalSecondDerivative(int part_idx, int direction, double h);
        /// \see Wavefunction::DriftForce(int part_idx)
        vector<double> DriftForce(int part_idx);
};