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

    // things I need multiple times in the calculations: better calculate them once for all
    double alpha = this->system->getWavefunction()->getParameter(0);
    double mass = 1.0; // !!!!!!!!!!!!!!!! Hard-code this is true only for the chosen model
    double res = 0.0;
    double psi = this->system->getWavefunction()->evaluate();
    double psidotdot = this->system->getWavefunction()->evaluateSecondDerivative();

    // potential contribution
    for(int i=0; i<this->system->getNParticles() ; i++){    
        for(int j=0; j<this->system->getDimension(); j++){
            res = res + pow( this->system->particles[i]->position[j], 2);
        }
    } 

    res = res * (0.5 * pow(this->omega,2) * mass - 2*pow(alpha, 2) / mass);
    res = res + this->system->getNParticles() * this->system->getDimension() * alpha/mass;

    double EL2 = 1/psi * ((-0.5)*psidotdot + this->potential()*psi);

    cout << res << " " << EL2 << endl;

    return res;

    // second derivative contribution
    //res = res - (0.5/this->system->particles[0]->mass/this->system->getWavefunction()->evaluate())*this->system->getWavefunction()->evaluateSecondDerivative();
}
