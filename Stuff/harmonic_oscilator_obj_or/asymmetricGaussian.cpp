#include "asymmetricGaussian.h"

AsymmetricGaussian::AsymmetricGaussian(System* s, double alpha, double beta) : Wavefunction(s){
    assert(s->getDimension() == 3); // Works only in 3D
    this->alpha = alpha;
    this->beta = beta;
}

double AsymmetricGaussian::evaluate(){
    double sum = 0.0; // Keeps track of the sum inside the exponetial
    /* Important!!!!
    When we have a sum all the constants that are repeatedly summed (e.g. alpha in this case) should be left outside of the sum
    and added in the multiplication at the end. In this wave we make the code faster because the require only one step to be added
    in the sum (otherwise N). */
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<2; j++){
            sum += pow(this->s->particles[i]->getPosition().at(j), 2);
        }
        sum += beta * pow(this->s->particles[i]->getPosition().at(2), 2); // Asymmetric therm
    }
    sum *= -this->alpha;
    return exp(sum);
}

double AsymmetricGaussian::evaluateSecondDerivative(){
    double wavef = evaluate();
    return (4*pow(this->alpha, 2) - 2*this->alpha) * wavef;
}

double AsymmetricGaussian::numericalSecondDerivative() {return 0.0;}