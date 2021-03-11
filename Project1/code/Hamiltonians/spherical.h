#pragma once
#include "hamiltonian.h"
#include <cmath>

class Spherical: public Hamiltonian{
    public:

        // Constructor and destructor
        Spherical(class System* system, double omega);
        ~Spherical();

        // Other functions
        /// \see LocalEnergyNumeric(double h), Hamiltonian::LocalEnergyAnalytic()
        double LocalEnergyAnalytic();
        /// \see LocalEnergyAnalytic(), Hamiltonian::LocalEnergyNumeric()
        double LocalEnergyNumeric(double h);
};
