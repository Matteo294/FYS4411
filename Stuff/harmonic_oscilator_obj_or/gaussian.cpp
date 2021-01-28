#include "gaussian.h"

Gaussian::Gaussian(System* s, double alpha) : Wavefunction(s, alpha){
    this->alpha = alpha;
}

double Gaussian::evaluate(){
    double sum = 0.0; // Keeps track of the sum inside the exponetial
    /* Important!!!!
    When we have a sum all the constants that are repeatedly summed (e.g. alpha in this case) should be left outside of the sum
    and added in the multiplication at the end. In this wave we make the code faster because the require only one step to be added
    in the sum (otherwise N). */
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<this->s->getDimension(); j++){
            sum += pow(this->s->particles[i]->position[j], 2);
        }
    }
    sum *= -this->alpha;
    return exp(sum);
}

double Gaussian::evaluateSecondDerivative(){
    double res = 0.0;
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<this->s->getDimension(); j++){
            res = res + pow(this->s->particles[i]->getPosition().at(j), 2);
        }
    }
    res = res*4*pow(this->alpha, 2);
    for(int i=0; i<(this->s->getDimension()*this->s->getNParticles()); i++){
        res = res - 2*this->alpha;
    }
    return res*evaluate();
}

double Gaussian::numericalSecondDerivative() {return 0.0;}