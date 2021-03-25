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

        /** Initializes the system with a certain dimensionality and a certain number of particles. If the simulations are run
         * in parallel, bool parallel is true. Particles are initialized in the origin, they will be moved while calling the termalizer.
         * \see Solver::thermalize() **/
        System(int dim, int Npart, bool parallel);
        ~System();           

        /** Returns the pointer to the selected Hamiltonian **/
        class Hamiltonian* getHamiltonian();
        /** Returns the pointer to the selected Wavefunction **/
        class Wavefunction* getWavefunction();
        /** Returns the pointer to the selected Solver **/
        class Solver* getSolver();
        /** Returns the pointer to the random generator **/
        class RandomGenerator* getRandomGenerator();
        /** Returns the dimensionality of the system **/
        int getDimension();
        /** Returns the number of particles present in the system **/
        int getNParticles();
        /** Returns True if the evaluation of the matrices for relative distance and position is enabled for the current simulation.**/
        bool getUseMatrix();
        /** Returns True if the simulation is run in parallel. **/
        bool getParallel();
        /** Vector of pointers to Particles object 
         * \see class Particle) **/
        vector<class Particle*> getParticles();

        /** Sets the Hamiltonian for the system, \see class Hamiltonian **/
        void setHamiltonian(class Hamiltonian* hamiltonian);
        /** Sets the Solver for the system, \see class Solver **/
        void setSolver(class Solver* solver);
        /** Sets the Wavefunvtion for the system, \see class Wavefunction **/
        void setWavefunction(class Wavefunction* wavefunction);
        /** Sets the Random generator for the system, \see class RandomGenerator **/
        void setRandomGenerator(class RandomGenerator* randomgenerator);
        /** If set to True, imposes the update of the relative_positions and relative_distances matrices every time that a particle is moved
         * This flag is automatically assigned when we set the wavefunction to the system and the choice is made depending on the number
         * of paramerers characterizing the assigned wavefunction:
         * Gaussian does not require the mentioned matrices and has 1 parameter, while AsymmGaussian gaussian requires it and has
         * 2 parameters. 
         * \see System::setWavefunction() **/
        void setUseMatrix(bool usematrix);

        /** Adds a particle to the system, placing it in the position specified by pos.**/
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
        class Hamiltonian* hamiltonian;
        class Wavefunction* wavefunction;
        class Solver* solver;
        vector<class Particle*> particles;  
        class RandomGenerator* randomgenerator;
        int dimension;
        int Nparticles;
        // Flag to activate/deactivate the relative positions and relative distances matrices
        bool usematrix;
        bool parallel;
        
};
