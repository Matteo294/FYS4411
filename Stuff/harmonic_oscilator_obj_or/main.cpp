#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "system.h"
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
    const int dimension = 1;
    const int Nparticles = 1;

    System system(dimension, Nparticles);
    Gaussian psi(&system, alpha);
    system.setWavefunction(&psi);

    //Spherical spherical(&system, omega);
    //system.setHamiltonian(&spherical);

    // generic info on the syste,
    system.getInfoParticles();
    cout << system.getNParticles() << " " << system.getDimension() << endl;

    // evaluate the wf with all the particles in the origin
    cout << "gaussian wf = " <<  system.getWavefunction()->evaluate() << endl;
    // change the position of two particles
    vector<double> new_pos {1}; // Chosen so that r1=r2=1 -> expected psi(r1=1, r2=2) = 1/e^2 = 0.135
    system.setParticlePosition(0, new_pos);
    // re-evaluate the wavefunction
    cout << "gaussian wf = " <<  system.getWavefunction()->evaluate() << endl;

    // local energy for the system
    cout << "Localenergy = " << system.hamiltonian->LocalEnergy() << endl;
}