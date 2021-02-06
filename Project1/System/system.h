#pragma once 
#include <vector>
#include "../Particles/particle.h"
#include "../Wavefunctions/wavefunction.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Solvers/solver.h"
#include "../Others/random_generator.h"
#include <iostream>
#include <cassert> 


using namespace std;

class System{
    public:

        // Constructor
        System(int dim, int Npart);             

        // Getters
        class Hamiltonian* getHamiltonian();
        class Wavefunction* getWavefunction();
        class Solver* getSolver();
        class RandomGenerator* getRandomGenerator();
        int getDimension();
        int getNParticles();
        vector<class Particle*> getParticles();


        // Setters
        void setHamiltonian(class Hamiltonian* hamiltonian);
        void setSolver(class Solver* solver);
        void setWavefunction(class Wavefunction* wavefunction);
        void setRandomGenerator(class RandomGenerator* randomgenerator);

        // Other functions
        void addParticle(double mass, vector<double> pos);

    private:
        // Attributes
        class Hamiltonian* hamiltonian;
        class Wavefunction* wavefunction;
        class Solver* solver;
        vector<class Particle*> particles; // vector of pointers to the single particles in the system   
        class RandomGenerator* randomgenerator;

        // Access this data only via getters
        int dimension; 
        int Nparticles;
};
