#include "asymmetricGaussian.h"

AsymmetricGaussian::~AsymmetricGaussian(){};

AsymmetricGaussian::AsymmetricGaussian(System* s, double alpha, double beta, double a) : Wavefunction(s, 3){
    assert(s->getDimension() == 3); // Works only in 3D
    this->setParameter(0, alpha);
    this->setParameter(1, beta);
    this->setParameter(2, a);
    this->s->setUseMatrix(true);
}

double AsymmetricGaussian::evaluateAll(){
    double arg = 0.0; // Keeps track of the sum inside the exponetial
    double f = 1.0; // this accounts for the interaction term that multipies the gaussian
    double a = this->params[2];
    double relative_dist=0.0;
    int k=0, j=0;
    
    for(int i=0; i<this->s->getNParticles(); i++){
        for(k=0; k<i; k++){
            if(this->s->relative_distance[i][k]<a){
                f=0.0; 
            }
            f *= (1 - this->params[2] / this->s->relative_distance[i][k]);
        }
        
        for(j=0; j<2; j++){
            arg += (double) pow(this->s->getParticles()[i]->getPosition().at(j), 2);
        }
        arg += (double) this->getParameter(1) * pow(this->s->getParticles()[i]->getPosition().at(2), 2); // Asymmetric therm
    }
    arg *= -this->getParameter(0);
    return exp(arg) * f;
}

double AsymmetricGaussian::evaluateSing(int part_idx){
    double arg = 0.0;
    double f = 1.0;
    double relative_dist = 0.0;
    double a = this->s->getWavefunction()->getParameter(2);

    for(int i=0; i<this->s->getNParticles(); i++){
        
        if(i!=part_idx){
            if(this->s->relative_distance[part_idx][i]<a){
                f=0.0;
            }
            f *= 1 - this->params[2] / this->s->relative_distance[part_idx][i];
        }
    }

    for(int i=0; i<2; i++){
        arg += pow( this->s->getParticles()[part_idx]->getPosition()[i], 2);
    }
    arg += this->getParameter(1) * pow(this->s->getParticles()[part_idx]->getPosition()[2] , 2);
    return exp( -this->getParameter(0) * arg ) * f;
}

/* need to be fixed */
vector<double> AsymmetricGaussian::DriftForce(int part_idx) {
    vector<double> x{3,0}; 
    return x;
}

/* need to be fixed */
double AsymmetricGaussian::analyticalAlphaDerivative(){
    double res=0.0;
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<2; j++){
            res += pow(this->s->getParticles()[i]->getPosition()[j], 2);
        }
        res += this->getParameter(1)*pow(this->s->getParticles()[i]->getPosition()[2], 2);
    }
    return -res * this->evaluateAll();
}

/* need to be fixed */
double AsymmetricGaussian::numericalSecondDerivative(int part_idx, int direction, double h) {return 0.0;}