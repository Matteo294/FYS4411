#pragma once 
#include "wavefunction.h"
#include <cmath>

class Gaussian: public Wavefunction{
    public:
        Gaussian(class System* s, double alpha);
        ~Gaussian();

        /// \see Wavefunction::evaluateAll()
        double evaluateAll();
        /// \see Wavefunction::evaluateSing(int part_idx)
        double evaluateSing(int part_idx);
        /// \see Wavefunction::numericalSecondDerivative(int part_idx, int direction, double h)
        double numericalSecondDerivative(int part_idx, int direction, double h);
        /// \see Wavefunction::analyticalAlphaDerivative()
        double analyticalAlphaDerivative(); 
        /// \see Wavefunction::DriftForce(int part_idx)
        vector<double> DriftForce(int part_idx);
};