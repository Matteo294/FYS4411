#include "gaussian.h"

Gaussian::Gaussian(System* s, double alpha) : Wavefunction(s){
    this->alpha = alpha;
}

double Gaussian::evaluate(){
    double sum = 0.0; // Keeps track of the sum inside the exponetial
    /* Important!!!!
    When we have a sum all the constants that are repeatedly summed (e.g. alpha in this case) should be left outside of the sum
    and added in the multiplication at the end. In this wave we make the code faster because the require only one step to be added
    in the sum (otherwise N). */
    for(int i=0; i<s->Nparticles; i++){
        for(int j=0; j<s->dimension; j++){
            sum += pow(s->particles[i]->position[j], 2);
        }
    }
    sum *= -alpha;
    return exp(sum);
}

double Gaussian::evaluateSecondDerivative(){
    double wavef = evaluate();
    return (4*pow(alpha, 2) - 2*alpha) * wavef;
}

double Gaussian::numericalSecondDerivative() {return 0.0;}