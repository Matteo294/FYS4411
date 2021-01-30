#pragma once
#include "system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *system);
        class System* system;

        /* uses analytic expression of the second derivative of the wavefunction to evaluate the result.
        At the moment this is implemented only for the non-interacting case */
        virtual double LocalEnergySecondDerivative() = 0;

        /* uses analytic expression of the local energy to evaluate the result.
        At the moment this is implemented only for the non-interacting case */
        virtual double LocalEnergyAnalytic() = 0;

};
