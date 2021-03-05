#include "spherical.h"

Spherical::Spherical(System* system, double omega) : Hamiltonian(system){
    this->omega=omega;
}


double Spherical::LocalEnergyAnalytic(){
    // things I need multiple times in the calculations or long expressions: better calculate them once for all
    // !!!! Hard-code! this is true only for the chosen model
    double alpha = this->system->getWavefunction()->getParameter(0);
    double mass = this->system->getParticles()[0]->getMass(); 
    int i=0, j=0;
    double res = 0;

    res = this->system->r2((double) 1.0);
    res *= 0.5 * mass * pow(this->omega, 2) - 2 * pow(alpha, 2) / mass;

    return res + alpha * this->system->getDimension() * this->system->getNParticles() / mass;

}



double Spherical::LocalEnergyNumeric(double h){
    // things I need multiple times in the calculations or long expressions: better calculate them once for all
    double alpha = this->system->getWavefunction()->getParameter(0);
    double mass = this->system->getParticles()[0]->getMass(); // !!!!!!!!!!!!!!!! Hard-code this is true only for the chosen model
    int i=0, j=0;
    double res = 0.0;
    
    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            res += this->system->getWavefunction()->numericalSecondDerivative(i, j, h);
        }
    }

    res *= -0.5 / mass / this->system->getWavefunction()->evaluateAll();
    
    res += 0.5 * mass * pow(this->omega, 2) * this->system->r2((double) 1.0);
    return res;
}
