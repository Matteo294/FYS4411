#include "particle.h"
#include "system.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include "wavefunction.h"
#include "gaussian.h"

using namespace std;

double rando(double min, double max){
    return min + ((double) rand()/RAND_MAX)*(max-min);
}

int main(){
    double min = -1.0;
    double max = 1.0;
    double alpha = 1.0;

    class System system;
    int dimension = 3;
    int Nparticles = 2;
    system.Initialize(dimension, Nparticles);
    system.getInfoParticles();
    cout << system.getNParticles() << " " << system.getDimension() << endl;
    /*system.addParticle();
    system.particles.back()->setPosition({1,2,3});
    system.particles.back()->setMass(1);
    system.getInfoParticles();*/

    vector<double> new_pos {sqrt(0.2), sqrt(0.5), sqrt(0.3)}; // Chosen so that r1=r2=1 -> expected psi(r1=1, r2=2) = 1/e^2 = 0.135
    system.particles[0]->setPosition(new_pos);
    system.particles[1]->setPosition(new_pos);
    Gaussian psi(&system, alpha);
    cout << "Value of psi: " << psi.evaluate() << endl;    

}