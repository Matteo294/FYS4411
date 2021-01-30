#include "metropolis.h"

Metropolis::Metropolis(System* system, int Nsteps, double step, double initialFraction) : Solver(system, Nsteps, step, initialFraction){

}

vector<double> Metropolis::solve(){

    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;
    double old_val = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = 2*this->step*( this->system->randomuniform() - 0.5);
        }

        this->system->setParticlePosition(i, new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->randomuniform() * (this->system->getNParticles() - 1));

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = 2*this->step*( this->system->randomuniform() - 0.5);
        }
        
        old_val = this->system->getWavefunction()->evaluateSing(idx);
        this->system->moveParticle(idx, new_pos);

        if( this->system->randomuniform() > (pow(this->system->getWavefunction()->evaluateSing(idx), 2) / pow(old_val, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->moveParticle(idx, {new_pos});
        }

        if(i>=(int)(this->Nsteps*this->initialFraction)){
            tmp = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
            energy += tmp;
            energy2 += tmp*tmp;
        }
        
    }

    //cout << "Energy = " << energy/this->Nsteps/(1-this->initialFraction) << endl;
    //cout << "STD = " << energy2/this->Nsteps/(1-this->initialFraction) - energy*energy/this->Nsteps/(1-this->initialFraction)/this->Nsteps/(1-this->initialFraction) << endl;
    return {energy/this->Nsteps/(1-this->initialFraction), energy2/this->Nsteps/(1-this->initialFraction) - energy*energy/this->Nsteps/(1-this->initialFraction)/this->Nsteps/(1-this->initialFraction)};
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