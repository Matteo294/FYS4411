#include "system.h"

/* ASSERT: Here I check that the dimension and the # of particles and if they satisfy the requirement
then I proceed with the initialization.
I initialize all the particles in the origin (with proper dimensions), each mass is set to 1 */

System::System(int dim, int Npart) : relative_position(Npart), relative_distance(Npart){
    assert(dim>0 && Npart>=0);
    int i=0;
    this->dimension=dim;
    this->Nparticles=Npart;
    this->usematrix = true; // this will be adjusted within the wavefunction class

    vector<double> pos(dim, 0.0);
    for(int i=0; i<Npart; i++){
        this->relative_distance[i].resize(this->Nparticles, 0);
        this->relative_position[i].resize(this->Nparticles, pos);
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
    res += parameter * pow(vect[this->dimension - 1], 2);
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


// Matrices
void System::EvaluateRelativePosition(){

    vector<double> rel_pos(this->dimension, 0.0);
    vector<double> pos_i(this->dimension, 0.0);
    vector<double> pos_j(this->dimension, 0.0);         

    for(int i=0; i<this->Nparticles; i++){
        pos_i = this->getParticles()[i]->getPosition();
        
        for(int j=0; j<i; j++){
            pos_j = this->getParticles()[j]->getPosition();
            transform(pos_i.begin(), pos_i.end(), pos_j.begin(), rel_pos.begin(), minus<double>());
            this->relative_position[i][j] = rel_pos;
            transform(rel_pos.begin(), rel_pos.end(), rel_pos.begin(), bind1st(multiplies<double>(), -1.0));
            this->relative_position[j][i] = rel_pos;
        }

        fill(this->relative_position[i][i].begin(), this->relative_position[i][i].end(), 0.0);
    }

}

void System::EvaluateRelativePosition(int idx){
    vector<double> rel_pos(this->dimension, 0.0);
    vector<double> pos_i(this->dimension, 0.0);
    vector<double> pos_j(this->dimension, 0.0);         
    pos_i = this->getParticles()[idx]->getPosition();

    for(int j=0; j<this->Nparticles; j++){
        if(j==idx){
            fill(this->relative_position[j][j].begin(), this->relative_position[j][j].end(), 0.0);
        } else {
            pos_j = this->getParticles()[j]->getPosition();
            transform(pos_i.begin(), pos_i.end(), pos_j.begin(), rel_pos.begin(), minus<double>());
            this->relative_position[idx][j] = rel_pos;
            transform(rel_pos.begin(), rel_pos.end(), rel_pos.begin(), bind1st(multiplies<double>(), -1.0));
            this->relative_position[j][idx] = rel_pos;
        }
    }

}

void System::EvaluateRelativeDistance(){
    for(int i=0; i<Nparticles; i++){
        for(int j=0; j<i; j++){
            this->relative_distance[i][j] = sqrt(this->r2(this->relative_position[i][j], (double) 1.0));
            this->relative_distance[j][i] = this->relative_distance[i][j];
        }
    }
}


void System::EvaluateRelativeDistance(int idx){
    for(int i=0; i<Nparticles; i++){
        this->relative_distance[idx][i] = sqrt(this->r2(this->relative_position[idx][i], (double) 1.0));
        this->relative_distance[i][idx] = this->relative_distance[idx][i];
    }
}


// Getters
        class Hamiltonian* System::getHamiltonian() {return this->hamiltonian;}
        class Wavefunction* System::getWavefunction() {return this->wavefunction;}
        class Solver* System::getSolver(){return this->solver;}
        class RandomGenerator* System::getRandomGenerator(){ return this->randomgenerator; }
        int System::getDimension() {return this->dimension;}
        int System::getNParticles() {return this->Nparticles;}
        vector<class Particle*> System::getParticles(){return this->particles;}
        bool System::getUseMatrix(){ return this->usematrix;}
        
        

// Setters
        void System::setHamiltonian(class Hamiltonian* hamiltonian) {this->hamiltonian = hamiltonian;}
        void System::setSolver(class Solver* solver) {this->solver = solver;}
        void System::setWavefunction(class Wavefunction* wavefunction) {this->wavefunction = wavefunction;}
        void System::setRandomGenerator(class RandomGenerator* randomgenerator) { this->randomgenerator = randomgenerator; }
        void System::setUseMatrix(bool usematrix){ this->usematrix = usematrix; }