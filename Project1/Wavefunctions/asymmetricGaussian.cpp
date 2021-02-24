#include "asymmetricGaussian.h"

AsymmetricGaussian::AsymmetricGaussian(System* s, double alpha, double beta) : Wavefunction(s, 2){
    assert(s->getDimension() == 3); // Works only in 3D
    this->setParameter(0, alpha);
    this->setParameter(1, beta);
    this->alpha = alpha;
    this->beta = beta;
}

double AsymmetricGaussian::evaluateAll(){
    double arg = 0.0; // Keeps track of the sum inside the exponetial
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<2; j++){
            arg += pow(this->s->getParticles()[i]->getPosition().at(j), 2);
        }
        arg += beta * pow(this->s->getParticles()[i]->getPosition().at(2), 2); // Asymmetric therm
    }
    arg *= -this->alpha;
    return exp(arg);
}

double AsymmetricGaussian::evaluateSing(int part_idx){
    double arg = 0.0;
    for(int i=0; i<2; i++){
        arg += pow(this->s->getParticles()[part_idx]->getPosition()[i], 2);
    }
    arg += this->beta * pow( this->s->getParticles()[part_idx]->getPosition()[2], 2);
    arg *= -this->alpha;
    return exp(arg);
}

/* !!!!!!!!!!!!! wrong, needs to be fixed !!!!!!!!!!!!!!!!!!! */
/*double AsymmetricGaussian::analyticalSecondDerivative(){
    double wavef = evaluateAll();
    double sum = 0.0;
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<2; j++){
            sum += (double) (4*pow(this->alpha, 2)*pow(this->s->getParticles()[i]->getPosition()[j], 2) - 2*this->alpha);
        }
        sum += (double) (4*pow(this->alpha*this->beta, 2)*pow(this->s->getParticles()[i]->getPosition()[2], 2) - 2*this->alpha*this->beta);
    }
    return sum*wavef;
}*/

vector<double> AsymmetricGaussian::DriftForce(int part_idx) {
    vector<double> x{3,0}; 
    return x;
}

double AsymmetricGaussian::analyticalAlphaDerivative(){
    double res=0.0;
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<2; j++){
            res += pow(this->s->getParticles()[i]->getPosition()[j], 2);
        }
        res += this->beta*pow(this->s->getParticles()[i]->getPosition()[2], 2);
    }
    return -res * this->evaluateAll();
}


double AsymmetricGaussian::numericalSecondDerivative(int part_idx, int direction, double h) {return 0.0;}