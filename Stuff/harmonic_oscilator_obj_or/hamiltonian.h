#pragma once
#include "system.h"

class Hamiltonian{
    public:
        virtual double LocalEnergy()=0;
        Hamiltonian(class System *system);

    protected:
        class System* system;
};
