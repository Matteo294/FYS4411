#include "spherical.h"

Spherical::Spherical(System* system, double omega) : Hamiltonian(system){
    this->omega=omega;
}

double Spherical::potential(){
    double x2 = 0.0;
    for(int i=0; i<this->system->getNParticles(); i++){
        for(int j=0; j<this->system->getDimension(); j++){
            x2 += pow(this->system->particles[i]->getPosition().at(j), 2); // xi^2 + yi^2 + zi^2 for each particle
        }
    }
    return 0.5 * pow(this->omega, 2) * x2;
}

double Spherical::LocalEnergy(){

    // things I need multiple times in the calculations or long expressions: better calculate them once for all
    double alpha = this->system->getWavefunction()->getParameter(0);
    double mass = 1.0; // !!!!!!!!!!!!!!!! Hard-code this is true only for the chosen model
    double res = 0.0;
    double psi = this->system->getWavefunction()->evaluate();
    double psidotdot = this->system->getWavefunction()->evaluateSecondDerivative();

    // Local energy
    double EL = 1/psi * ((-0.5)*psidotdot + this->potential()*psi);

    return EL;
}
