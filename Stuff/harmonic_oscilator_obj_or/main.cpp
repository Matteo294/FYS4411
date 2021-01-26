#include "particle.h"
#include "system.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <cstdlib>
using namespace std;

double rando(double min, double max){
    return min + ((double) rand()/RAND_MAX)*(max-min);
}

int main(){
    double min = -1.0;
    double max = 1.0;

    class System system;
    int dimension = 3;
    int Nparticles = 100;
    system.Initialize(dimension, Nparticles);
    system.getInfoParticles();
    cout << system.getNParticles() << " " << system.getDimension() << endl;
    system.addParticle();
    system.particles.back()->setPosition({1,2,3});
    system.particles.back()->setMass(1);
    system.getInfoParticles();
    

}