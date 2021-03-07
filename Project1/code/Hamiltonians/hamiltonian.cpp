#include "hamiltonian.h"

Hamiltonian::~Hamiltonian(){};

Hamiltonian::Hamiltonian(System* system){
    this->system = system;
}
