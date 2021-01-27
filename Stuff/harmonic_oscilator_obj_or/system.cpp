#include "system.h"

/* ASSERT: Here I check that the dimension and the # of particles and if they satisfy the requirement
then I proceed with the initialization.
I initialize all the particles in the origin (with proper dimensions), each mass is set to 1 */

System::System(int dim, int Npart){
    assert(dim>0 && Npart>=0);
    int i=0, j=0;
    this->dimension=dim;
    this->Nparticles=Npart;

    vector<double> pos;
    for(i=0; i<this->dimension; i++){
        pos.push_back((double) 0.0);
    }

    for(i=0; i<Npart; i++){
        this->addParticle();
        for(j=0; j<this->dimension; j++){
            particles[i]->setMass((double) 1.0);
            particles[i]->setPosition(pos);
        }
    }
}

void System::setHamiltonian(class Hamiltonian* hamiltonian){
    this->hamiltonian=hamiltonian;
}

class Hamiltonian* System::getHamiltonian(){
    return this->hamiltonian;
}

void System::setWavefunction(class Wavefunction* wavefunction){
    this->wafefunction=wavefunction;
}

class Wavefunction* System::getWavefunction(){
    return this->wafefunction;
}

/* Add a particle by adding a new pointer to the vector of pointer */
void System::addParticle(){
    this->particles.push_back(new Particle());
    this->Nparticles=particles.size();    
}

int System::getDimension(){
    return this->dimension; 
}

int System::getNParticles(){
    return this->Nparticles;
}

void System::getInfoParticles(){
    int i=0;
    int j=0;
    for(i=0; i<this->Nparticles; i++){
        cout << "Particle[" << i << "]   mass=" << particles[i]->mass << "   Position={";
        for(j=0; j<this->dimension; j++){
            cout << this->particles[i]->getPosition().at(j);
            if(j != (this->dimension -1)){
                cout << ",";
            }
        }
        cout << "}" << endl;
    }
}

void System::moveParticle(int i, vector<double> delta_pos){
    particles[i]->move(delta_pos);
}

void System::setParticlePosition(int i, vector<double> new_pos){
    assert(i < Nparticles);
    particles[i]->setPosition(new_pos);
}