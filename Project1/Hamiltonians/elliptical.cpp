#include "elliptical.h"

Elliptical::Elliptical(System* system, double omegaXY, double omegaZ) : Hamiltonian(system){
    this->omegaXY = omegaXY;
    this->omegaZ = omegaZ;
}

double Elliptical::potential(){
    // Hardcoded mass=1
    double V = 0.0;
    for(int i=0; i<this->system->getNParticles(); i++){
        for(int j=0; j<2; j++){
            V += pow(this->getOmegaXY()*this->system->getParticles()[i]->getPosition()[j], 2);
        }
        V += pow(this->getOmegaZ()*this->system->getParticles()[i]->getPosition()[2], 2);
    }
    return 0.5*V;
}

double Elliptical::LocalEnergyAnalytic(){
    // things I need multiple times in the calculations or long expressions: better calculate them once for all
    // Hardcoded mass=1
    double alpha = this->system->getWavefunction()->getParameter(0);
    int i=0, j=0;
    double res = 0;

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<2; j++){
            res += pow(this->system->getParticles()[i]->getPosition()[j], 2);
        }
        res += pow(this->system->getWavefunction()->getParameter(1) * this->system->getParticles()[i]->getPosition()[2], 2);        
    }

    return 3*alpha*this->system->getNParticles() - 2*pow(alpha, 2)*res + this->potential();

}

double Elliptical::LocalEnergyNumeric(double h){return 0.0;}

/*double Elliptical::LocalEnergyNumeric(double h){
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
    
    double pot = 0.0;
    for(int i=0; i<this->system->getNParticles(); i++){
        for(int j=0; j<this->system->getDimension(); j++){
            pot += pow(this->system->getParticles()[i]->getPosition().at(j), 2); // xi^2 + yi^2 + zi^2 for each particle
        }
    }
    res += 0.5 * mass * pow(this->omega, 2) * pot;
    return res;
}*/
