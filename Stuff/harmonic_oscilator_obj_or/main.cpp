#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "system.h"
#include "gaussian.h"
#include "asymmetricGaussian.h"
#include "spherical.h"
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
    const int Nparticles = 3;

    System system(dimension, Nparticles);
    Gaussian psi(&system, alpha);
    Spherical spherical(&system, omega);

    system.setWavefunction(&psi);
    system.setHamiltonian(&spherical);

    // generic info on the syste,
    system.getInfoParticles();
    cout << system.getNParticles() << " " << system.getDimension() << endl;
    vector<double> new_pos {1,0,0}; // Chosen so that r1=r2=1 -> expected psi(r1=1, r2=2) = 1/e^2 = 0.135
    system.setParticlePosition(0, new_pos);
    system.getInfoParticles();
    cout << "arrivo qui" << endl;
    
    cout << "Local Energy = " << system.getHamiltonian()->LocalEnergy() << endl;
    cout << "gaussian wf = " <<  system.getWavefunction()->evaluate() << endl;
}