#include "hamiltonian.h"

Hamiltonian::~Hamiltonian(){};

Hamiltonian::Hamiltonian(System* system, int nparams){
    this->system = system;
    this->nparams = nparams;
    this->params.resize(nparams, 0.0);
}


void Hamiltonian::setParameter(int idx, double value){
    this->params[idx] = value;
}


int Hamiltonian::getnparameter(){ return this->nparams; }
double Hamiltonian::getParameter(int idx){
    assert(idx < this->nparams);
    return this->params[idx];
}
