#include "wavefunction.h"

Wavefunction::~Wavefunction(){};

Wavefunction::Wavefunction(System *s, int nparams){
    this->s = s; 
    this->params.resize(nparams, 0.0); 
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

int Wavefunction::getNparams(){ return this->nparams; }
