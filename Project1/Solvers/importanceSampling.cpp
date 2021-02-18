#include "importanceSampling.h"
#include <fstream>

ImportanceSampling::ImportanceSampling(System* system, int Nsteps, double initialFraction, double dt, double D) : Solver(system, Nsteps, initialFraction){
    this->dt=dt;
    this->D = D;
}

vector<double> ImportanceSampling::computeEnergy(){
    random_device rd;
    mt19937_64 gen(rd());
        
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;

    double old_val = 0.0;
    double arg = 0.0;
    
    vector<double> new_pos(this->system->getDimension(), 0.0);
    vector<double> old_drift(this->system->getDimension(), 0.0);
    vector<double> new_drift(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->dt);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        old_val = this->system->getWavefunction()->evaluateSing(idx);
        old_drift = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->D * old_drift[j] * this->dt + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }
        
        double arg = 0.0;
        this->system->getParticles()[idx]->move(new_pos);
        new_drift = this->system->getWavefunction()->DriftForce(idx);
        for(j=0; j<this->system->getDimension(); j++){
            arg += (old_drift[j] + new_drift[j]) * ( -2 * new_pos[j] + this->D * this->dt * (old_drift[j] - new_drift[j]) );
        }

        arg *= 0.25;
        

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(this->system->getWavefunction()->evaluateSing(idx), 2) / pow(old_val, 2))){
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


vector<double> ImportanceSampling::computeEnergy(double h){
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;

    double old_val = 0.0;
    double arg = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);
    vector<double> old_drift(this->system->getDimension(), 0.0);
    vector<double> new_drift(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->dt);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        old_val = this->system->getWavefunction()->evaluateSing(idx);
        old_drift = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->D * old_drift[j] * this->dt + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }
        
        this->system->getParticles()[idx]->move(new_pos);
        new_drift = this->system->getWavefunction()->DriftForce(idx);
        for(j=0; j<this->system->getDimension(); j++){
            arg += (old_drift[j] + new_drift[j]) * ( -2 * new_pos[j] + this->D * this->dt * (old_drift[j] - new_drift[j]) );
        }

        arg *= 0.25;
        

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(this->system->getWavefunction()->evaluateSing(idx), 2) / pow(old_val, 2))){
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

vector<double> ImportanceSampling::computeExpectationPsidotPsi(){
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double x=0.0, x2=0.0, tmp=0.0;

    double old_val=0.0, new_val=0.0;
    double psidot=0.0;
    double arg = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);
    vector<double> old_drift(this->system->getDimension(), 0.0);
    vector<double> new_drift(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->dt);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        old_val = this->system->getWavefunction()->evaluateSing(idx);
        old_drift = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->D * old_drift[j] * this->dt + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }
        
        this->system->getParticles()[idx]->move(new_pos);

        // Updates drift force and wavefunction
        new_drift = this->system->getWavefunction()->DriftForce(idx);
        new_val = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            arg += (old_drift[j] + new_drift[j]) * ( -2 * new_pos[j] + this->D * this->dt * (old_drift[j] - new_drift[j]) );
        }

        arg *= 0.25;
        
        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(new_val, 2) / pow(old_val, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->getParticles()[idx]->move(new_pos);
        } else{
            new_val = old_val; // useful for the following lines
        }

        psidot = this->system->getWavefunction()->analyticalAlphaDerivative();
        if(i>=(int)(this->Nsteps*this->InitialFraction)){
            tmp = (double) psidot / new_val;
            x += tmp;
            x2 += tmp*tmp;
        }
        
    }

    return {x/this->Nsteps/(1-this->InitialFraction), x2/this->Nsteps/(1-this->InitialFraction) - x*x/this->Nsteps/(1-this->InitialFraction)/this->Nsteps/(1-this->InitialFraction)};
}

vector<double> ImportanceSampling::computeExpectationPsidotPsiEl(){
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double x=0.0, x2=0.0, tmp=0.0;

    double old_val=0.0, new_val=0.0;
    double psidot=0.0;
    double arg = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);
    vector<double> old_drift(this->system->getDimension(), 0.0);
    vector<double> new_drift(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->dt);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        old_val = this->system->getWavefunction()->evaluateSing(idx);
        old_drift = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->D * old_drift[j] * this->dt + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }
        
        this->system->getParticles()[idx]->move(new_pos);

        // Updates drift force and wavefunction
        new_drift = this->system->getWavefunction()->DriftForce(idx);
        new_val = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            arg += (old_drift[j] + new_drift[j]) * ( -2 * new_pos[j] + this->D * this->dt * (old_drift[j] - new_drift[j]) );
        }

        arg *= 0.25;
        
        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(new_val, 2) / pow(old_val, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->getParticles()[idx]->move(new_pos);
        } else{
            new_val = old_val; // useful for the following lines
        }

        psidot = this->system->getWavefunction()->analyticalAlphaDerivative();
        if(i>=(int)(this->Nsteps*this->InitialFraction)){
            tmp = (double) psidot / new_val * system->getHamiltonian()->LocalEnergyAnalytic();
            x += tmp;
            x2 += tmp*tmp;
        }
        
    }

    return {x/this->Nsteps/(1-this->InitialFraction), x2/this->Nsteps/(1-this->InitialFraction) - x*x/this->Nsteps/(1-this->InitialFraction)/this->Nsteps/(1-this->InitialFraction)};
}


// Gettes
    double ImportanceSampling::getdt(){ return this->dt; }
    double ImportanceSampling::getD(){ return this->D; }

// Setters
    void ImportanceSampling::setdt(double dt) { this->dt = dt; }
    void ImportanceSampling::setD(double D) { this->D = D; }
