#pragma once 
#include <vector>
#include "particle.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "solver.h"
#include <iostream>
#include <cassert> 

using namespace std;

class System{
    public:
        class Hamiltonian* hamiltonian;
        void setHamiltonian(class Hamiltonian* hamiltonian);
        class Hamiltonian* getHamiltonian();

        class Wavefunction* wafefunction;
        void setWavefunction(class Wavefunction* wavefunction);
        class Wavefunction* getWavefunction();

        class Solver* solver;
        void setSolver(class Solver* solver);
        class Solver* getSolver();

        vector<class Particle*> particles; // vector of pointers to the single particles in the system
        void addParticle(double mass, vector<double> pos);
        int getDimension();
        int getNParticles();
        void getInfoParticles();
        void moveParticle(int i, vector<double> delta_pos);
        void setParticlePosition(int i, vector<double> new_pos);

        System(int dim, int Npart);
    private:
        int dimension; 
        int Nparticles;
};
