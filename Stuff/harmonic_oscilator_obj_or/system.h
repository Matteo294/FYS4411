#pragma once 
#include <vector>
#include "particle.h"
#include "wavefunction.h"
#include "hamiltonian.h"
#include "solver.h"
#include "random_generator.h"
#include <iostream>
#include <cassert> 


using namespace std;

class System{
    public:

        // Constructor
        System(int dim, int Npart);

        // Attributes
        class Hamiltonian* hamiltonian;
        class Wavefunction* wavefunction;
        class Solver* solver;
        vector<class Particle*> particles; // vector of pointers to the single particles in the system   
        class RandomGenerator* randomgenerator;
             

        // Getters
        class Hamiltonian* getHamiltonian();
        class Wavefunction* getWavefunction();
        class Solver* getSolver();
        class RandomGenerator* getRandomGenerator();
        int getDimension();
        int getNParticles();
        void getInfoParticles();


        // Setters
        void setParticlePosition(int i, vector<double> new_pos);
        void setHamiltonian(class Hamiltonian* hamiltonian);
        void setSolver(class Solver* solver);
        void setWavefunction(class Wavefunction* wavefunction);
        void setRandomGenerator(class RandomGenerator* randomgenerator);

        // Other functions
        void addParticle(double mass, vector<double> pos);
        void moveParticle(int i, vector<double> delta_pos);

    private:
        // Access this data only via getters
        int dimension; 
        int Nparticles;
};
