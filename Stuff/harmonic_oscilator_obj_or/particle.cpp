#include "particle.h"

Particle::Particle(System* system, double mass, std::vector<double> pos){
    this->position = pos;
    this->mass = mass;
    this->system=system;
}

void Particle::setMass(double m){
    this->mass = m;
}

void Particle::setPosition(std::vector<double> new_pos){
    this->position = new_pos;
}

void Particle::move(std::vector<double> var_pos){
    assert(var_pos.size()==this->system->getDimension());
    for(int i=0; i<this->system->getDimension(); i++){
        this->position.at(i) += var_pos.at(i);
    }
}

std::vector<double> Particle::getPosition(){
    return this->position;
}