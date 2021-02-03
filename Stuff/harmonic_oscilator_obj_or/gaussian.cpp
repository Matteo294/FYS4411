#include "gaussian.h"
#include <vector>

Gaussian::Gaussian(System* s, double alpha) : Wavefunction(s, 1){
    this->setParameter(0, alpha);
    this->alpha = alpha; // <-- Is this necessary?
}

double Gaussian::evaluateAll(){
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

double Gaussian::evaluateSing(int part_idx){
    double arg = 0.0;
    for(int i=0; i<this->s->getDimension(); i++){
        arg += pow( this->s->particles[part_idx]->getPosition()[i] , 2);
    }
    return exp( -this->alpha * arg );
}

double Gaussian::evaluateSecondDerivative(){
    double res = 0.0;
    for(int i=0; i<this->s->getNParticles(); i++){
        for(int j=0; j<this->s->getDimension(); j++){
            res = res + pow(this->s->particles[i]->getPosition().at(j), 2); // particle i at its position 3 element vector alla dimensione j  
        }
    }
    res = res*4*pow(this->alpha, 2);
    for(int i=0; i<(this->s->getDimension()*this->s->getNParticles()); i++){
        res = res - 2*this->alpha;
    }
    //std::cout << "res1 " << evaluateAll() << endl;
    return res*evaluateAll();
}

vector<double> Gaussian::evaluateGradient(){
    vector<double> grad;
    // mac isn't able to recognise this command grad.reshape(3,0.0);
    grad.push_back(0);
    grad.push_back(0);
    grad.push_back(0);
    double ev=evaluateAll();
    for(int j=0; j<this->s->getDimension(); j++){ // for cicle on all the dimentions
        for(int i=0; i<this->s->getNParticles(); i++){ // for cicle on all the particles
            grad[j] += -2*this->alpha*this->s->particles[i]->getPosition().at(j);   
        }
        grad[j] *= ev;
    }
    //std::cout << " Gradient:  " << endl << "x: " << grad[0] << " y: " << grad[1]/evaluateAll()  << " z: " << grad[2] << endl;
    //std::cout << " position0: " << endl << "x:" << this->s->particles[0]->getPosition().at(0) << " y:" << this->s->particles[0]->getPosition().at(1) << " z:" << this->s->particles[0]->getPosition().at(2)<< endl;
    //std::cout << " position1: " << endl << "x:" << this->s->particles[1]->getPosition().at(0) << " y:" << this->s->particles[1]->getPosition().at(1) << " z:" << this->s->particles[1]->getPosition().at(2)<< endl;
    return grad;
}


double Gaussian::numericaldriftForce(int part_idx, int direction, double h){ //input: particle ID, direction where calculate the derivative, increment 
    assert( direction < this->s->getDimension() );
    assert ( part_idx < this->s->getNParticles() );
    double res = 0.0;
    int i=0, j=0;

    vector<double> drift(this->s->getDimension(), 0.0);
    // evaluate wf(x+h)
    double wf_start =evaluateAll(); 
    drift[direction] = h;
    this->s->moveParticle(part_idx, drift);
    double wf_forw = evaluateAll();

    // final result
    drift[direction] = -h;
    this->s->moveParticle(part_idx, drift); // particles returns to its original position
    res =wf_forw - wf_start;
    //std::cout << "Gradient numerical in x: " << res/h << endl;
    return res/h;

}

double Gaussian::numericalSecondDerivative(int part_idx, int direction, double h) {
    assert( direction < this->s->getDimension() );
    assert ( part_idx < this->s->getNParticles() );
    double res = 0.0;
    int i=0, j=0;

    // initialize a vector for moving the particle
    vector<double> drift(this->s->getDimension(), 0.0);
    
    // evaluate wf(x+h)
    drift[direction] = h;
    this->s->moveParticle(part_idx, drift);
    double wf_forw = evaluateAll();

    // evaluate wf(x-h)
    drift[direction] = -2*h;
    this->s->moveParticle(part_idx, drift);
    double wf_back = evaluateAll();

    // final result
    drift[direction] = h;
    this->s->moveParticle(part_idx, drift); // particles returns to its original position
    res = wf_back + wf_forw - 2*evaluateAll();

    return res/pow(h, 2);

}