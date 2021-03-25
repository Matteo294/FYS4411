#include "gaussian.h"

Gaussian::~Gaussian(){};

Gaussian::Gaussian(System* s, double alpha) : Wavefunction(s, 1){
    this->setParameter(0, alpha);
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

    // Initialize a vector for moving the particle
    vector<double> pos_var(this->s->getDimension(), 0.0);
    
    // evaluate wf(x+h)
    pos_var[direction] = h;
    this->s->getParticles()[part_idx]->move(pos_var);
    double wf_forw = this->evaluateSing(part_idx);

    // evaluate wf(x-h)
    pos_var[direction] = -2*h;
    this->s->getParticles()[part_idx]->move(pos_var);
    double wf_back = this->evaluateSing(part_idx);

    // final result
    pos_var[direction] = h;
    this->s->getParticles()[part_idx]->move(pos_var); // particle returns to its original position
    res = wf_back + wf_forw - 2 * this->evaluateSing(part_idx);

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