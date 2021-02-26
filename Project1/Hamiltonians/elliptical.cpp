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
    double beta = this->system->getWavefunction()->getParameter(1);
    double a = this->system->getWavefunction()->getParameter(2);
    int i=0, j=0, m=0;
    double res = 0.0;
    vector<double> vec(this->system->getDimension(), 0.0);
    vector<double> vec2(this->system->getDimension(), 0.0);
    double tmp=0.0;

    for(i=0; i<this->system->getNParticles(); i++){
        
        for(m=0; m<this->system->getNParticles(); m++){
            if(m!=i){
                vec = this->system->getParticles()[i]->getRelativePosition(m);
                //cout << vec[0] << vec[1] << vec[2] << endl;
                tmp = this->system->getParticles()[i]->getRelativeDistance(m);
                //cout << tmp << endl;
                transform(vec.begin(), vec.end(), vec.begin(), bind1st(multiplies<double>(), a / tmp / (tmp - a)));
            }

            transform(vec2.begin(), vec2.end(), vec.begin(), vec2.begin(), plus<double>());
        }

    }

    cout << vec2[0] << " " << vec2[1] << " " << vec2[2] << endl;

    return 0.0;
}

double Elliptical::LocalEnergyNumeric(double h){return 0.0;}
