#pragma once
#include "hamiltonian.h"
#include <cmath>

class Spherical: public Hamiltonian{
    public:
        Spherical(class System* system, double omega);
        double omega;
        double LocalEnergy();
};
