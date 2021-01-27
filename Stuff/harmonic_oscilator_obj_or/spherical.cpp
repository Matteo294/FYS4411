#include "spherical.h"

Spherical::Spherical(System* system, double omega) : Hamiltonian(system){
    this->omega=omega;
}

double Spherical::LocalEnergy(){
    double res = 0.0;
    // potential contribution
    for(int i=0; i<this->system->getNParticles() ; i++){    
        for(int j=0; j<this->system->getDimension(); j++){
            res = res + pow( this->system->particles[i]->position[j], 2);
        }
    }     

    res = res*0.5*pow(this->omega,2)*this->system->particles[0]->mass;

    // second derivative contribution
    res = res - (0.5/this->system->particles[0]->mass/this->system->getWavefunction()->evaluate())*this->system->getWavefunction()->evaluateSecondDerivative();

    return res;
}

