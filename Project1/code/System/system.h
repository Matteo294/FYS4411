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
#include <omp.h>


using namespace std;

class System{
    public:

        // Constructor and destructor
        System(int dim, int Npart, bool parallel);
        ~System();           

        // Getters
        class Hamiltonian* getHamiltonian();
        class Wavefunction* getWavefunction();
        class Solver* getSolver();
        class RandomGenerator* getRandomGenerator();
        int getDimension();
        int getNParticles();
        vector<class Particle*> getParticles();
        bool getUseMatrix();
        bool getParallel();

        // Setters
        void setHamiltonian(class Hamiltonian* hamiltonian);
        void setSolver(class Solver* solver);
        void setWavefunction(class Wavefunction* wavefunction);
        void setRandomGenerator(class RandomGenerator* randomgenerator);
        /// \see relative_position, relative_distance
        void setUseMatrix(bool usematrix);

        // Other functions
        /// Adds a particle to the system
        void addParticle(double mass, vector<double> pos); 
        /** Evaluates the squared distance of particles and sums.
         * The parameter can be set to give a different coefficient to the last squared coordinate 
         * of each particle (useful for asymmetric wavefunction and elliptical potential)
         * \note set parameter to one if you want a spherical potential and a symmetrical wavefunction
        **/
        double r2(double parameter);
        /** Evaluates the sum of the square of the vector components.
         * The parameter gives a different coefficient to the last squared coordinate of vect
         * \note set parameter to one if you want to use the simple sum of the squared components
        **/
        double r2(vector<double> vect, double parameter);
        /// Scalar product between vector 1 and vector 2
        double cdot(vector<double> v1, vector<double> v2);

        // Matrices
        /// NxN matrix of 3d vectors. The ij-th element is a 3D vector containing posi - posj
        vector<vector<vector<double>>> relative_position;
        /// NxN matrix of doubles. The ij-th element is the distance between particle i and particle j
        vector<vector<double>> relative_distance;
        /// Updates the relative_distance matrix. \see relative_distance
        void EvaluateRelativeDistance();
        /// Updates the relative position matrix. \see relative_position
        void EvaluateRelativePosition();
        /// Updates a special row of the relative position matrix. \see relative_position
        void EvaluateRelativePosition(int idx);
        /// Updates a special row of the relative distance matrix. \see relative_distance
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
        // Flag to activate/deactivate the relative positions and relative distances matrices
        bool usematrix;
        bool parallel;
        
};
