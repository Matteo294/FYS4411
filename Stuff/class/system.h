#ifndef Sasha_H
#define Sasha_H

#include "particle.h"
#include <vector>

using namespace std;

class System{
    private:
        int n = 0; // keeps track of the number of particles
        vector<Particle> particles; // This allows us to create the space in the memory for the objects without running constructors (we cannot predict values :) )
    public:
        void addParticle(Particle p);   // This function should be used to add a Particle to the Sasha
        void accelerateAll(); // Accelerates all particles
        void getInfo(); // Info of all particles
};

#endif