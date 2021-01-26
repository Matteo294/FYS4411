#pragma once 
#include <vector>
#include "particle.h"
#include <iostream>
#include <cassert> // library for assert function (used in system.cpp)

using namespace std;

class System{
    public:
        vector<class Particle*> particles; // vector of pointers to the single particles in the system
        int dimension; 
        int Nparticles;
        void addParticle();
        int getDimension();
        int getNParticles();
        void getInfoParticles();
        void moveParticle(int i, vector<double> delta_pos);
        void setParticlePosition(int i, vector<double> new_pos);
        System(int dim, int Npart);
};