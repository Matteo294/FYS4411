#include "asymmetricGaussian.h"

AsymmetricGaussian::AsymmetricGaussian(System* s, double alpha, double beta) : Wavefunction(s, 2){
    assert(s->getDimension() == 3); // Works only in 3D
    this->setParameter(0, alpha);
    this->setParameter(1, beta);
    this->alpha = alpha;
    this->beta = beta;
}

double AsymmetricGaussian::evaluateAll(){
    double sum = 0.0; // Keeps track of the sum inside the exponetial
    /* Important!!!!
    When we have a sum all the constants that are repeatedly summed (e.g. alpha in this case) should be left outside of the sum
    and added in the multiplication at the end. In this wave we make the code faster because the require only one step to be added
    in the sum (otherwise N). */
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<2; j++){
            sum += pow(this->s->getParticles()[i]->getPosition().at(j), 2);
        }
        sum += beta * pow(this->s->getParticles()[i]->getPosition().at(2), 2); // Asymmetric therm
    }
    sum *= -this->alpha;
    return exp(sum);
}

double AsymmetricGaussian::evaluateSing(int part_idx){
    double arg = 0.0;
    for(int i=0; i<2; i++){
        arg += pow( this->s->getParticles()[part_idx]->getPosition()[i], 2);
    }
    arg += this->beta * pow( this->s->getParticles()[part_idx]->getPosition()[2], 2);
    arg *= -this->alpha;

    return exp(arg);
}

/* !!!!!!!!!!!!! wrong, needs to be fixed !!!!!!!!!!!!!!!!!!! */
double AsymmetricGaussian::evaluateSecondDerivative(){
    double wavef = evaluateAll();
    return (4*pow(this->alpha, 2) - 2*this->alpha) * wavef;
}

double AsymmetricGaussian::numericalSecondDerivative(int part_idx, int direction, double h) {return 0.0;}