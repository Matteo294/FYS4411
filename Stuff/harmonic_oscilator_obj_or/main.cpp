#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>

#include "particle.h"
#include "system.h"
#include "gaussian.h"
#include "asymmetricGaussian.h"

using namespace std;

double rando(double min, double max){
    return min + ((double) rand()/RAND_MAX)*(max-min);
}

int main(){

    const double min = -1.0;
    const double max = 1.0;
    const double alpha = 1.0;
    const double beta = 1.0;
    const int dimension = 3;
    const int Nparticles = 2;

    System system(dimension, Nparticles);
    
    system.getInfoParticles();
    cout << system.getNParticles() << " " << system.getDimension() << endl;

    vector<double> new_pos {sqrt(0.2), sqrt(0.5), sqrt(0.3)}; // Chosen so that r1=r2=1 -> expected psi(r1=1, r2=2) = 1/e^2 = 0.135
    system.setParticlePosition(0, new_pos);
    system.setParticlePosition(1, new_pos);

    // Compare results between symmetric gaussian and asymmetric gaussian but with alpha=beta
    Gaussian psi(&system, alpha);
    AsymmetricGaussian psi2(&system, alpha, beta);
    cout << "Value of psi: " << psi.evaluate() << endl;
    cout << "Valule of psi2: " << psi2.evaluate() << endl;

}