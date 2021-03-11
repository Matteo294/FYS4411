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
    int k=0, j=0;
    
    for(int i=0; i<this->s->getNParticles(); i++){
        for(k=0; k<i; k++){
            if(this->s->relative_distance[i][k]<this->params[2]){
                f=0.0; 
            }
            f *= (1 - this->params[2] / this->s->relative_distance[i][k]);
        }
    }
    
    return exp( - this->params[0] * this->s->r2((double) this->params[1])) * f;
}

double AsymmetricGaussian::evaluateSing(int part_idx){
    double arg = 0.0;
    double f = 1.0;
    double relative_dist = 0.0;

    for(int i=0; i<this->s->getNParticles(); i++){
        if(i!=part_idx){
            if(this->s->relative_distance[part_idx][i]<this->params[2]){
                f=0.0;
            }
            f *= 1 - this->params[2] / this->s->relative_distance[part_idx][i];
        }
    }
    
    return exp( - this->params[0] * this->s->r2(this->s->getParticles()[part_idx]->getPosition(), (double) this->params[1]) ) * f;
}

/* need to be fixed */
vector<double> AsymmetricGaussian::DriftForce(int part_idx) {
    // reducing the notation, we will call part_idx as i in nominanting vectors
    vector<double> r_ij(3, 0.0);
    vector<double> res(3, 0.0);
    double dist_ij=0.0;

    vector<double> r_i = this->s->getParticles()[part_idx]->getPosition();
    r_i[2] *= this->params[1]; // I multiply the third element by beta
    transform(r_i.begin(), r_i.end(), r_i.begin(), bind1st(multiplies<double>(), - 2 * this->params[0])); // I multiply the vector by -2*alpha
    
    for(int j=0; j<this->s->getNParticles(); j++){
        if(j!=part_idx){
            r_ij = this->s->relative_position[part_idx][j];
            dist_ij = this->s->relative_distance[part_idx][j];
            transform(r_ij.begin(), r_ij.end(), r_ij.begin(), bind1st(multiplies<double>(), this->params[2] / (dist_ij - this->params[2]) / pow(dist_ij, 2) ) );
            transform(res.begin(), res.end(), r_ij.begin(), res.begin(), plus<double>());
        }
    }
    //at the end of the cycle I evaluated the whole sum over j appearing in the analytical expression
    // now I add the last term
    transform(res.begin(), res.end(), r_i.begin(), res.begin(), plus<double>());
    transform(res.begin(), res.end(), res.begin(), bind1st(multiplies<double>(), 2.0) );
    return res;
}


double AsymmetricGaussian::psibar_psi(){
    return -this->s->r2(this->params[1]);
}


double AsymmetricGaussian::numericalSecondDerivative(int part_idx, int direction, double h) {
    assert( direction < this->s->getDimension() );
    assert ( part_idx < this->s->getNParticles() );
    double res = 0.0;

    // initialize a vector for moving the particle
    vector<double> pos_var(this->s->getDimension(), 0.0);
    
    // evaluate wf(x+h)
    pos_var[direction] = h;
    this->s->getParticles()[part_idx]->move(pos_var);
    this->s->EvaluateRelativePosition(part_idx); this->s->EvaluateRelativeDistance(part_idx);
    double wf_forw = this->evaluateAll();

    // evaluate wf(x-h)
    pos_var[direction] = -2*h;
    this->s->getParticles()[part_idx]->move(pos_var);
    this->s->EvaluateRelativePosition(part_idx); this->s->EvaluateRelativeDistance(part_idx);
    double wf_back = this->evaluateAll();

    // final result
    pos_var[direction] = h;
    this->s->getParticles()[part_idx]->move(pos_var); // particles returns to its original position
    this->s->EvaluateRelativePosition(part_idx); this->s->EvaluateRelativeDistance(part_idx);
    res = wf_back + wf_forw - 2 * this->evaluateAll();

    return res/pow(h, 2);
}