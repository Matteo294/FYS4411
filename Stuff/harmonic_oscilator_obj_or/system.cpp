#include "system.h"

System::System(){}

/* ASSERT: Here I check that the dimension and the # of particles and if they satisfy the requirement
then I proceed with the initialization.
I initialize all the particles in the origin (with proper dimensions), each mass is set to 1 */

void System::Initialize(int dim, int Npart){
    assert(dim>0 && Npart>=0);
    int i=0, j=0;
    this->dimension = dim;

    vector<double> pos;
    for(i=0; i<this->dimension; i++){
        pos.push_back((double) 0.0);
    }

    for(i=0; i<Npart; i++){
        this->addParticle();
        for(j=0; j<this->dimension; j++){
            this->particles.at(i)->setMass((double) 1.0);
            this->particles.at(i)->setPosition(pos);
        }
    }

    this->Nparticles=particles.size();

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
        cout << "Particle[" << i << "]   mass=" << this->particles.at(i)->mass << "   Position={";
        for(j=0; j<this->dimension; j++){
            cout << this->particles.at(i)->getPosition().at(j);
            if(j!= (this->dimension -1)){
                cout << ",";
            }
        }
        cout << "}" << endl;
    }
}