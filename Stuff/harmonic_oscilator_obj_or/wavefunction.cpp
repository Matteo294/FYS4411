#include "wavefunction.h"

Wavefunction::Wavefunction(System *s, double alpha){
    this->s = s; // Simply saves the system pointer in the class
    this->alpha=alpha;
}