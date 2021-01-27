#include "metropolis.h"

Metropolis::Metropolis(System* system, int Nsteps, double step, double initialFraction) : Solver(system, Nsteps, step, initialFraction){

}

void Metropolis::solve(){
    srand(100);
    int i=0, j=0, idx=0;
    double energy=0, energy2 = 0, tmp=0;
    vector<double> new_pos;

    //initialization in random positions
    for(i=0; i<system->getDimension(); i++){
        new_pos.push_back((double) 0.0);
    }

    for(i=0; i<system->getNParticles(); i++){
        for(j=0; j<system->getDimension(); j++){
            new_pos[j] = this->step*(((double) rand()/RAND_MAX) - 0.5);
        }

        system->setParticlePosition(i, new_pos);
    }

    
    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round(((double) rand()/RAND_MAX)*(system->getNParticles() - 1));

        for(j=0; j<system->getDimension(); j++){
            new_pos[j] = this->step*(((double) rand()/RAND_MAX) - 0.5);
        }
        
        system->setParticlePosition(idx, new_pos);

        /*
        if(i>(this->Nsteps*this->initialFraction)){
            tmp = system->getHamiltonian()->LocalEnergy();
            energy += tmp;
            energy2 += tmp*tmp;
        }
        */

        
        

    }

    
    cout << "Energy = " << energy/this->Nsteps << endl;
    cout << "STD = " << energy2/this->Nsteps - energy*energy/this->Nsteps/this->Nsteps << endl;
    
}

int Metropolis::getNsteps(){
    return this->Nsteps;
}

double Metropolis::getstep(){
    return this->step;
}

double Metropolis::getinitialFraction(){
    return this->initialFraction;
}