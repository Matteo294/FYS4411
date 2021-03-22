#pragma once
#include "hamiltonian.h"
#include <cmath>

class Spherical: public Hamiltonian{
    public:

        // Constructor and destructor
        Spherical(class System* system, double omega);
        ~Spherical();

        // Other functions
        double LocalEnergyAnalytic();
        double LocalEnergyNumeric(double h);
};
