#pragma once 
#include <vector>
#include "../Particles/particle.h"
#include "../Wavefunctions/wavefunction.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Solvers/solver.h"
#include "../Others/random_generator.h"
#include "../Others/functions.h"
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
        bool getUseMatrix();

        // Setters
        void setHamiltonian(class Hamiltonian* hamiltonian);
        void setSolver(class Solver* solver);
        void setWavefunction(class Wavefunction* wavefunction);
        void setRandomGenerator(class RandomGenerator* randomgenerator);
        void setUseMatrix(bool usematrix);

        // Other functions
        void addParticle(double mass, vector<double> pos); 
        double r2(double parameter);
        double r2(vector<double> vect, double parameter); 
        double cdot(vector<double> v1, vector<double> v2);

        // Matrices
        vector<vector<vector<double>>> relative_position;
        vector<vector<double>> relative_distance;
        void EvaluateRelativeDistance();
        void EvaluateRelativePosition();
        void EvaluateRelativePosition(int idx);
        void EvaluateRelativeDistance(int idx);
        
        

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
        bool usematrix;
        
};
