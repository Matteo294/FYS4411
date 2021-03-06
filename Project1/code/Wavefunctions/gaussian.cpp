#include "gaussian.h"

Gaussian::Gaussian(System* s, double alpha) : Wavefunction(s, 1){
    this->setParameter(0, alpha);
    this->s->setUseMatrix(false);
}

double Gaussian::evaluateAll(){ 
    return exp(-this->params[0] * this->s->r2((double) 1.0));
}

double Gaussian::evaluateSing(int part_idx){
    return exp(-this->params[0] * this->s->r2(this->s->getParticles()[part_idx]->getPosition(), (double) 1.0));
}

double Gaussian::numericalSecondDerivative(int part_idx, int direction, double h) {
    assert( direction < this->s->getDimension() );
    assert ( part_idx < this->s->getNParticles() );
    double res = 0.0;

    // initialize a vector for moving the particle
    vector<double> drift(this->s->getDimension(), 0.0);
    
    // evaluate wf(x+h)
    drift[direction] = h;
    this->s->getParticles()[part_idx]->move(drift);
    double wf_forw = evaluateAll();

    // evaluate wf(x-h)
    drift[direction] = -2*h;
    this->s->getParticles()[part_idx]->move(drift);
    double wf_back = evaluateAll();

    // final result
    drift[direction] = h;
    this->s->getParticles()[part_idx]->move(drift); // particles returns to its original position
    res = wf_back + wf_forw - 2*evaluateAll();

    return res/pow(h, 2);

}

vector<double> Gaussian::DriftForce(int part_idx){
    vector<double> v = this->s->getParticles()[part_idx]->getPosition();
    transform(v.begin(), v.end(), v.begin(), bind1st(multiplies<double>(), -4 * this->params[0]));  
    return v;
}


double Gaussian::psibar_psi(){
    return  - this->s->r2((double) 1.0);
}