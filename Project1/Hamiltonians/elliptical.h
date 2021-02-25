#pragma once
#include "hamiltonian.h"
#include <cmath>

class Elliptical: public Hamiltonian{
    public:

        // Constructor
        Elliptical(class System* system, double omegaXY, double omegaZ);

        // Getters
        double getOmegaXY() {return this->omegaXY;}
        double getOmegaZ() {return this->omegaZ;}

        //Other functions
        double LocalEnergyAnalytic();
        double LocalEnergyNumeric(double h);
        double potential();

    private:
        // Private attributes: access only via getters and setters
        double omegaXY;
        double omegaZ;
};
