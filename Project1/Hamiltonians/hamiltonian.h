#pragma once
#include "../System/system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *system);
        class System* system;

        /* uses analytic expression of the local energy to evaluate the result.
        At the moment this is implemented only for the non-interacting case */
        virtual double LocalEnergyAnalytic() = 0;
        virtual double LocalEnergyNumeric(double h) = 0;
};
