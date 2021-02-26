#include "system.h"

/* ASSERT: Here I check that the dimension and the # of particles and if they satisfy the requirement
then I proceed with the initialization.
I initialize all the particles in the origin (with proper dimensions), each mass is set to 1 */

System::System(int dim, int Npart){
    assert(dim>0 && Npart>=0);
    int i=0;
    this->dimension=dim;
    this->Nparticles=Npart;

    // we create all the particles in the origin
    vector<double> pos(dim, 0.0);

    for(i=0; i<this->dimension; i++){
        pos.push_back((double) 0.0);
    }

    this->particles.resize(Nparticles);

    for(int i=0; i<this->Nparticles; i++){
        this->particles[i] = new Particle(this, 1.0, pos);
    }
}

/* Add a particle by adding a new pointer to the vector of pointer - not tested yet */
void System::addParticle(double mass, vector<double> pos){
    this->particles.push_back(new Particle(this, mass, pos));
    this->Nparticles++;
}


/* returns (x^2 + y^2 + param * z^2) for particle part_idx*/
double System::r2(vector<double> vect, double parameter){
    double res = 0.0;
    for(int i=0; i<(this->dimension - 1); i++){
        res += pow(vect[i], 2);
    }
    res += parameter * pow(vect.back(), 2);

    return res;
}

/* returns the sum over all particles of (x^2 + y^2 + param * z^2) */
double System::r2(double parameter){
    double res = 0.0;
    vector<double> pos(this->dimension, 0.0);

    for(int i=0; i<this->Nparticles; i++){
        pos = this->getParticles()[i]->getPosition();
        res += this->r2(pos, parameter);
    }
    
    return res;
}


double System::cdot(vector<double> v1, vector<double> v2){
    double res = 0.0;
    for(int i=0; i<this->dimension; i++){
        res += v1[i] * v2[i];
    }

    return res;
}


// Getters
        class Hamiltonian* System::getHamiltonian() {return this->hamiltonian;}
        class Wavefunction* System::getWavefunction() {return this->wavefunction;}
        class Solver* System::getSolver(){return this->solver;}
        class RandomGenerator* System::getRandomGenerator(){ return this->randomgenerator; }
        int System::getDimension() {return this->dimension;}
        int System::getNParticles() {return this->Nparticles;}
        vector<class Particle*> System::getParticles(){return this->particles;}
        

// Setters
        void System::setHamiltonian(class Hamiltonian* hamiltonian) {this->hamiltonian = hamiltonian;}
        void System::setSolver(class Solver* solver) {this->solver = solver;}
        void System::setWavefunction(class Wavefunction* wavefunction) {this->wavefunction = wavefunction;}
        void System::setRandomGenerator(class RandomGenerator* randomgenerator) { this->randomgenerator = randomgenerator; }