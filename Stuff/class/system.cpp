#include <iostream>
#include "particle.h"
#include "system.h"

/* Adds particle to the System */
void System::addParticle(Particle p){
    this->particles.push_back(p); // push back puts the element inside the brackets at the end of the vector
    this->n++; // Increments particles counter
}

void System::getInfo(){
    for(int i=0; i<this->n; i++){
        this->particles[i].getInfo();
    }
}

void System::accelerateAll(){
    for(int i=0; i<this->n; i++){
        this->particles[i].accelerate();
    }
}