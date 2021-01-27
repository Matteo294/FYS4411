#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "system.h"
#include "gaussian.h"
#include "asymmetricGaussian.h"
#include "spherical.h"
#include "metropolis.h"
using namespace std;

double rando(double min, double max){
    return min + ((double) rand()/RAND_MAX)*(max-min);
}

int main(){

    const double min = -1.0;
    const double max = 1.0;
    const double alpha = 1.0;
    const double beta = 1.0;
    double omega = 1.0;
    const int dimension = 3;
    const int Nparticles = 10;

    const int Nsteps = (int) 1e6;
    const double step = 1.0;
    const double initialFraction = 0.1;

    System system(dimension, Nparticles);
    Gaussian psi(&system, alpha);
    Spherical spherical(&system, omega);
    Metropolis metropolis(&system, Nsteps, step, initialFraction);

    system.setWavefunction(&psi);
    system.setHamiltonian(&spherical);
    system.setSolver(&metropolis);

    system.solver->solve();
    system.getInfoParticles();
    

}