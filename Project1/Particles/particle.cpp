#include "particle.h"


Particle::Particle(System* system, double mass, vector<double> pos){
    this->position = pos;
    this->mass = mass;
    this->system=system;
}

void Particle::setMass(double m){
    this->mass = m;
}

void Particle::setPosition(vector<double> new_pos){
    this->position = new_pos;
}

void Particle::move(vector<double> var_pos){
    assert(var_pos.size()==this->system->getDimension());
    for(int i=0; i<this->system->getDimension(); i++){
        this->position.at(i) += var_pos.at(i);
    }
}

vector<double> Particle::getPosition(){
    return this->position;
}

vector<double> Particle::getRelativePosition(int idx){
    vector<double> pos = this->system->getParticles()[idx]->getPosition(); // position of the other particle
    vector<double> res(this->system->getDimension(), 0.0); 
    
    transform(this->position.begin(), this->position.end(), pos.begin(), res.begin(), minus<double>());

    return res;
}

double Particle::getRelativeDistance(int idx){
    vector<double> rel_pos = this->getRelativePosition(idx);
    return sqrt(this->system->r2(rel_pos, (double) 1.0));
}