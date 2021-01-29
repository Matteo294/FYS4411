#pragma once
#include "hamiltonian.h"
#include <cmath>

class Spherical: public Hamiltonian{
    public:

        // Constructor
        Spherical(class System* system, double omega);

        // Getters
        double getOmega() {return this->omega;}

        //Other functions
        double LocalEnergy();
        double potential();

    private:
        // Private attributes: access only via getters and setters
        double omega;
};
