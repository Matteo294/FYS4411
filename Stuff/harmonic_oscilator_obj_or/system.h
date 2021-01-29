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

        System(int dim, int Npart); // Constructor

        // Attributes
        class Hamiltonian* hamiltonian;
        class Wavefunction* wavefunction;
        class Solver* solver;
        vector<class Particle*> particles; // vector of pointers to the single particles in the system        

        // Getters
        class Hamiltonian* getHamiltonian() {return this->hamiltonian;}
        class Wavefunction* getWavefunction() {return this->wavefunction;}
        class Solver* getSolver(){return this->solver;}
        int getDimension() {return this->dimension;}
        int getNParticles() {return this->Nparticles;}
        void getInfoParticles();

        // Setters
        void setParticlePosition(int i, vector<double> new_pos) {this->particles[i]->setPosition(new_pos);}
        void setHamiltonian(class Hamiltonian* hamiltonian) {this->hamiltonian = hamiltonian;}
        void setSolver(class Solver* solver) {this->solver = solver;}
        void setWavefunction(class Wavefunction* wavefunction) {this->wavefunction = wavefunction;}

        // Other functions
        void addParticle(double mass, vector<double> pos);
        void moveParticle(int i, vector<double> delta_pos) {this->particles[i]->move(delta_pos);}

    private:
        // Access this data only via getters
        int dimension; 
        int Nparticles;
};
