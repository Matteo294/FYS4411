#pragma once
#include "hamiltonian.h"
#include <cmath>

class Elliptical: public Hamiltonian{
    public:

        // Constructor and destructor
        Elliptical(class System* system, double omegaXY, double omegaZ);
        ~Elliptical();

        //Other functions
        double LocalEnergyAnalytic();
        double LocalEnergyNumeric(double h);

};
