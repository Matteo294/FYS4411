#include "metropolis.h"

Metropolis::Metropolis(System* system, int Nsteps, double initialFraction, double step) : Solver(system, Nsteps, initialFraction){
    this->step = step;
}

vector<double> Metropolis::solve(bool allAverages){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, psi_bar_psi=0.0, psi_bar_psi_EL=0.0, tmp1=0.0, tmp2=0.0;
    double psi_old = 0.0, psi_new=0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    bool last_accepted;
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        this->system->getParticles()[i]->setPosition(pos_var);
    }
    
    if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}


    //MCsteps
    //#pragma omp parallel for reduction(+:energy) reduction(+:energy2) \
    //firstprivate(gen, idx, pos_old, psi_old, pos_var, j, psi_new, last_accepted, tmp1, tmp2, psi_bar_psi, psi_bar_psi_EL) \
    //shared(accepted, allAverages)
    for(i=1; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));
        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }

        this->system->getParticles()[idx]->move(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
        psi_new = this->system->getWavefunction()->evaluateSing(idx);

        if( this->system->getRandomGenerator()->uniform(gen) > ( pow(psi_new,2) / pow(psi_old,2) )){
            this->system->getParticles()[idx]->setPosition(pos_old);
            if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
            last_accepted = 0;
        } else {
            accepted++;
            last_accepted = 1;
        }

        if( i>=(int)(this->Nsteps*this->InitialFraction) ){
            
            if(i==(int)(this->Nsteps*this->InitialFraction)){
                tmp1 = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
            } else {
                if(last_accepted){
                    tmp1 = this->system->getHamiltonian()->LocalEnergyAnalytic();
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
vector<double> Metropolis::solve(double h){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;
    double psi_old = 0.0, psi_new=0.0;
    int accepted = 0;
    bool usematrix = this->system->getUseMatrix();
    double ratio_accepted=0.0;
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    for(i=1; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        this->system->getParticles()[i]->setPosition(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->step*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        this->system->getParticles()[idx]->move(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
        psi_new = this->system->getWavefunction()->evaluateSing(idx);

        if( this->system->getRandomGenerator()->uniform(gen) > (pow(psi_new, 2) / pow(psi_old, 2))){
            this->system->getParticles()[idx]->setPosition(pos_old);
            if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
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

double Metropolis::getStep(){
    return this->step;
}
