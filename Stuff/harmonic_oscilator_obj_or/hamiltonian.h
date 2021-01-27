#pragma once
#include "system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *system);
        class System* system;
        virtual double LocalEnergy()=0;
        

};
