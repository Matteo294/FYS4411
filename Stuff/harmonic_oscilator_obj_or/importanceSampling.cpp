#include "importanceSampling.h"

ImportanceSampling::ImportanceSampling(System* system, int Nsteps, double step, double initialFraction, double dt, double D) : Solver(system, Nsteps, step, initialFraction){
    this->dt=dt;
    this->D = D;
}

vector<double> ImportanceSampling::solve(){
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
            new_pos[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }

        this->system->setParticlePosition(i, new_pos);
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        old_val = this->system->getWavefunction()->evaluateSing(idx);
        old_drift = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->getD() * old_drift[j] * this->getdt() + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }
        
        this->system->moveParticle(idx, new_pos);
        new_drift = this->system->getWavefunction()->DriftForce(idx);
        for(j=0; j<this->system->getDimension(); j++){
            arg += (old_drift[j] + new_drift[j]) * ( -2 * new_pos[j] + this->getD() * this->getdt() * (old_drift[j] - new_drift[j]) );
        }

        arg *= 0.25;
        

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(this->system->getWavefunction()->evaluateSing(idx), 2) / pow(old_val, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->moveParticle(idx, {new_pos});
        }

        if(i>=(int)(this->Nsteps*this->initialFraction)){
            tmp = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
            //tmp = (double) this->system->getHamiltonian()->LocalEnergyNumeric();
            energy += tmp;
            energy2 += tmp*tmp;
        }
        
    }

    //cout << "Energy = " << energy/this->Nsteps/(1-this->initialFraction) << endl;
    //cout << "STD = " << energy2/this->Nsteps/(1-this->initialFraction) - energy*energy/this->Nsteps/(1-this->initialFraction)/this->Nsteps/(1-this->initialFraction) << endl;
    return {energy/this->Nsteps/(1-this->initialFraction), energy2/this->Nsteps/(1-this->initialFraction) - energy*energy/this->Nsteps/(1-this->initialFraction)/this->Nsteps/(1-this->initialFraction)};
}

int ImportanceSampling::getNsteps(){
    return this->Nsteps;
}

double ImportanceSampling::getstep(){
    return this->step;
}

double ImportanceSampling::getinitialFraction(){
    return this->initialFraction;
}

double ImportanceSampling::getdt(){
    return this->dt;
}

double ImportanceSampling::getD(){
    return this->D;
}