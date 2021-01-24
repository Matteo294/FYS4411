#include "gaussian.h"
#include <cmath>

Gaussian::Gaussian(System* s, double alpha) : Wavefunction(s){
    this->alpha = alpha;
}

double Gaussian::evaluate(System* s){
    double sum = 0.0; // Keeps track of the sum inside the exponetial
    /* Important!!!!
    When we have a sum all the constants that are repeatedly summed (e.g. alpha in this case) should be left outside of the sum
    and added in the multiplication at the end. In this wave we make the code faster because the require only one step to be added
    in the sum (otherwise N). */
    for(int i=0; i < s->nparticles; i++){
        sum += pow(s.particles[i].x, 2);
    }
    sum *= -alpha;
    return exp(sum);
}

double Gaussian::evaluateSecondDerivative(System* s){
    double sum = 0.0;
    for(int i=0; i < s->nparticles; i++){
        sum += pow(s.particles[i].x, 2);
    }
    sum *= -alpha;
    return exp(sum);
}