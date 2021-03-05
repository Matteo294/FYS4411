#include "wavefunction.h"

Wavefunction::Wavefunction(System *s, int nparams){
    this->s = s; // Simply saves the system pointer in the class
    this->params.resize(nparams, 0.0); // Set number of parameters, set all parameters to 0
    this->nparams = nparams;
    
}

double Wavefunction::getParameter(int idx){ 
    assert(idx < this->nparams);
    return this->params[idx];
}

void Wavefunction::setParameter(int idx, double val){
    assert(idx < this->nparams);
    this->params[idx] = val;
}
