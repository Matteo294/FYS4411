#include "spherical.h"

Spherical::~Spherical(){};

Spherical::Spherical(System* system, double omega) : Hamiltonian(system, (int) 1){
    this->setParameter(0, omega);
}


double Spherical::LocalEnergyAnalytic(){
    // things I need multiple times in the calculations or long expressions: better calculate them once for all
    // !!!! Hard-code! this is true only for the chosen model
    double alpha = this->system->getWavefunction()->getParameter(0);
    double mass = this->system->getParticles()[0]->getMass(); 
    int i=0, j=0;
    double res = 0;

    res = this->system->r2((double) 1.0);
    res *= 0.5 * mass * pow(this->params[0], 2) - 2 * pow(alpha, 2) / mass;

    return res + alpha * this->system->getDimension() * this->system->getNParticles() / mass;

}



double Spherical::LocalEnergyNumeric(double h){
    double mass = this->system->getParticles()[0]->getMass(); 
    int i=0, j=0;
    double res_i = 0.0, res=0.0;
    
    for(i=0; i<this->system->getNParticles(); i++){
        res_i=0.0;
        for(j=0; j<this->system->getDimension(); j++){
            res_i += this->system->getWavefunction()->numericalSecondDerivative(i, j, h);
        }
        res_i /= this->system->getWavefunction()->evaluateSing(i);
        res += res_i;
    }

    res *= -0.5 / mass;
    res += 0.5 * mass * pow(this->params[0], 2) * this->system->r2((double) 1.0);
    return res;
    
}
