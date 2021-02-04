#include "metropolis.h"

Metropolis::Metropolis(System* system, int Nsteps, double initialFraction, double step) : Solver(system, Nsteps, initialFraction){
    this->step = step;
}

vector<double> Metropolis::solve(){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;
    double old_val = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        
        old_val = this->system->getWavefunction()->evaluateSing(idx);
        this->system->getParticles()[idx]->move(new_pos);

        if( this->system->getRandomGenerator()->uniform(gen) > (pow(this->system->getWavefunction()->evaluateSing(idx), 2) / pow(old_val, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->getParticles()[idx]->move(new_pos);
        }

        if(i>=(int)(this->Nsteps*this->InitialFraction)){
            tmp = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();            
            energy += tmp;
            energy2 += tmp*tmp;
        }
        
    }

    return {energy/this->Nsteps/(1-this->InitialFraction), energy2/this->Nsteps/(1-this->InitialFraction) - energy*energy/this->Nsteps/(1-this->InitialFraction)/this->Nsteps/(1-this->InitialFraction)};
}


vector<double> Metropolis::solve(double h){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;
    double old_val = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        
        old_val = this->system->getWavefunction()->evaluateSing(idx);
        this->system->getParticles()[idx]->move(new_pos);

        if( this->system->getRandomGenerator()->uniform(gen) > (pow(this->system->getWavefunction()->evaluateSing(idx), 2) / pow(old_val, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->getParticles()[idx]->move(new_pos);
        }

        if(i>=(int)(this->Nsteps*this->InitialFraction)){
            tmp = (double) this->system->getHamiltonian()->LocalEnergyNumeric(h);            
            energy += tmp;
            energy2 += tmp*tmp;
        }
        
    }

    return {energy/this->Nsteps/(1-this->InitialFraction), energy2/this->Nsteps/(1-this->InitialFraction) - energy*energy/this->Nsteps/(1-this->InitialFraction)/this->Nsteps/(1-this->InitialFraction)};
}


double Metropolis::getStep(){
    return this->step;
}
