#pragma once
#include "hamiltonian.h"
#include <cmath>

class Spherical: public Hamiltonian{
    public:

        // Constructor and destructor
        Spherical(class System* system, double omega);
        ~Spherical();

        // Getters
        double getOmega() {return this->omega;}

        // Other functions
        /// \see LocalEnergyNumeric(double h), Hamiltonian::LocalEnergyAnalytic()
        double LocalEnergyAnalytic();
        /// \see LocalEnergyAnalytic(), Hamiltonian::LocalEnergyNumeric()
        double LocalEnergyNumeric(double h);

    private:
        // Private attributes: access only via getters and setters
        double omega;
};
