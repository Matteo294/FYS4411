#include "particle.h"

Particle::Particle(){
}

void Particle::setMass(double m){
    this->mass = m;
}

void Particle::setPosition(std::vector<double> new_pos){
    this->position = new_pos;
}

void Particle::updPosition(std::vector<double> var_pos){
    assert(var_pos.size()==(*system).dimension);
    int i=0;
    for(i=0; i<(*system).dimension; i++){
        this->position.at(i) += var_pos.at(i);
    }
}

std::vector<double> Particle::getPosition(){
    return this->position;
}