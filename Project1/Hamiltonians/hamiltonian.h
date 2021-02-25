#pragma once
#include "../System/system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *system);
        class System* system;

        /* uses analytic expression of the local energy to evaluate the result.
        At the moment this is implemented only for the non-interacting case */
        virtual double LocalEnergyAnalytic() = 0;
        // returns the variation in the local energy after that the idx-th particle has been moved from old_pos to new_pos
        virtual double LocalEnergyVariation(int idx, vector<double> old_pos, vector<double> new_pos)=0;
        virtual double LocalEnergyNumeric(double h) = 0;

};
