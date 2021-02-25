#include "importanceSampling.h"
#include <fstream>

ImportanceSampling::ImportanceSampling(System* system, int Nsteps, double initialFraction, double dt, double D) : Solver(system, Nsteps, initialFraction){
    this->dt=dt;
    this->D = D;
}

vector<double> ImportanceSampling::solve(bool allAverages){
    random_device rd;
    mt19937_64 gen(rd());
        
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, psi_bar_psi=0.0, psi_bar_psi_EL=0.0, tmp1=0.0, tmp2=0.0;
    double psi_old = 0.0, psi_new=0.0;
    double arg = 0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    bool last_accepted;
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_new(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    vector<double> drift_old(this->system->getDimension(), 0.0);
    vector<double> drift_new(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->dt);
        }

        this->system->getParticles()[i]->setPosition(pos_var);
    }

    //MCsteps
    for(i=1; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->D * drift_old[j] * this->dt + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }

        this->system->getParticles()[idx]->move(pos_var);
        pos_new = this->system->getParticles()[idx]->getPosition();
        psi_new = this->system->getWavefunction()->evaluateSing(idx);
        drift_new = this->system->getWavefunction()->DriftForce(idx);
        
        arg = 0.0;
        for(j=0; j<this->system->getDimension(); j++){
            arg += (drift_old[j] + drift_new[j]) * ( -2 * pos_var[j] + this->D * this->dt * (drift_old[j] - drift_new[j]) );
        }

        arg *= 0.25;

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(psi_new, 2) / pow(psi_old, 2))){
            this->system->getParticles()[idx]->setPosition(pos_old);
            last_accepted = 0;
        } else {
            accepted++;
            last_accepted = 1;
        }

        if(i>=(int)(this->Nsteps*this->InitialFraction)){

             if(i==(int)(this->Nsteps*this->InitialFraction)){
                tmp1 = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
            } else {
                if(last_accepted){
                    tmp1 += this->system->getHamiltonian()->LocalEnergyVariation(idx, pos_old, pos_new);
                }
            }

            energy += tmp1;
            energy2 += tmp1*tmp1;

            if (allAverages){
                tmp2 = this->system->getWavefunction()->analyticalAlphaDerivative() / this->system->getWavefunction()->evaluateAll();
                psi_bar_psi += tmp2;
                psi_bar_psi_EL += tmp2 * tmp1;
            }
        }
    }

    energy = energy/this->Nsteps/(1-this->InitialFraction);
    energy2 = energy2/this->Nsteps/(1-this->InitialFraction);
    psi_bar_psi = psi_bar_psi/this->Nsteps/(1-this->InitialFraction);
    psi_bar_psi_EL = psi_bar_psi_EL/this->Nsteps/(1-this->InitialFraction);
    ratio_accepted = (double) accepted/this->Nsteps;

    return {energy, energy2 - pow(energy, 2), ratio_accepted, 2 * (psi_bar_psi_EL - psi_bar_psi * energy)};
}


/* here there are a LOT of things to fix, but we need to ask to morten before proceeding*/
vector<double> ImportanceSampling::solve(double h){
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    double psi_old = 0.0, psi_new = 0.0;
    double arg = 0.0;
    vector<double> new_pos(this->system->getDimension(), 0.0);
    vector<double> drift_old(this->system->getDimension(), 0.0);
    vector<double> drift_new(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->dt);
        }

        this->system->getParticles()[i]->setPosition(new_pos);
    }

    //MCsteps
    for(i=1; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            new_pos[j] = this->D * drift_old[j] * this->dt + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getdt());
        }
        
        this->system->getParticles()[idx]->move(new_pos);
        psi_new = this->system->getWavefunction()->evaluateSing(idx);
        drift_new = this->system->getWavefunction()->DriftForce(idx);
        
        for(j=0; j<this->system->getDimension(); j++){
            arg += (drift_old[j] + drift_new[j]) * ( -2 * new_pos[j] + this->D * this->dt * (drift_old[j] - drift_new[j]) );
        }

        arg *= 0.25;
        

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(psi_new, 2) / pow(psi_old, 2))){
            for(j=0; j<this->system->getDimension(); j++){
                new_pos[j] = -new_pos[j];
            }
            this->system->getParticles()[idx]->move(new_pos);
        } else {
            accepted++;
        }

        if(i>=(int)(this->Nsteps*this->InitialFraction)){
            tmp = (double) this->system->getHamiltonian()->LocalEnergyNumeric(h);
            energy += tmp;
            energy2 += tmp*tmp;
        }
        
    }

    energy = energy/this->Nsteps/(1-this->InitialFraction);
    energy2 = energy2/this->Nsteps/(1-this->InitialFraction);
    ratio_accepted = (double) accepted/this->Nsteps;

    return {energy, energy2 - pow(energy, 2), ratio_accepted};
}


// Gettes
    double ImportanceSampling::getdt(){ return this->dt; }
    double ImportanceSampling::getD(){ return this->D; }

// Setters
    void ImportanceSampling::setdt(double dt) { this->dt = dt; }
    void ImportanceSampling::setD(double D) { this->D = D; }
