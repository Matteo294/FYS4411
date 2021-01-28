#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>

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
    double alpha = 0.5;
    const double beta = 1.0;
    double omega = 1.0;
    const int dimension = 3;
    const int Nparticles = 100;

    const int Nsteps = (int) 1e4;
    const double step = 1.0;
    const double initialFraction = 0.1;

    System system(dimension, Nparticles);
    Gaussian psi(&system, alpha);
    Spherical spherical(&system, omega);
    Metropolis metropolis(&system, Nsteps, step, initialFraction);

    system.setWavefunction(&psi);
    system.setHamiltonian(&spherical);
    system.setSolver(&metropolis);
    
    
    double alpha_min = 0.5;
    double alpha_max = 5;
    int alpha_N = 1;
    ofstream fout("results.csv");
    fout << "alpha,energy,STD" << endl;
    for(int i=0; i<alpha_N; i++){
        alpha = alpha_min + i*(alpha_max - alpha_min)/alpha_N;
        system.setWavefunction(new Gaussian(&system, alpha));
        fout << alpha << "," << system.solver->solve()[0] << "," << system.solver->solve()[1] << endl;
    }

    fout.close();
    
    

}