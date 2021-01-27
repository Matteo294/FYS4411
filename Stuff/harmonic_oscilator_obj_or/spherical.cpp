#include "spherical.h"

Spherical::Spherical(System* system, double omega) : Hamiltonian(system){
    this->omega=omega;
}

double Spherical::LocalEnergy(){
    double res = 0.0;

    // potential contribution
    for(int i=0; i<system->getNParticles(); i++){    
        for(int j=0; j<system->getDimension(); j++){
            res = res + pow(system->particles[i]->position[j], 2);
        }
    }     
    res = res*0.5*omega*2*system->particles[0]->mass;

    // second derivative contribution
    res = res + (0.5/system->particles[0]->mass)*system->wafefunction->evaluateSecondDerivative();
    // final result
    res = system->wafefunction->evaluate()*res;

    return res;
}

