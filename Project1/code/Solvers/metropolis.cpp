#include "metropolis.h"
#include <omp.h>

Metropolis::~Metropolis(){};

Metropolis::Metropolis(System* system, int Nsteps, double initialFraction, double step) : Solver(system, Nsteps, initialFraction, 1){
    this->setParameter(0, step);
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
    int thermalization_step = (int) this->Nsteps*this->InitialFraction;
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    
    /* I DON'T KNOW IF THIS IS WORKING PROPERLY
    #pragma omp parallel for collapse(2)
    for(j=0; j<this->system->getDimension(); j++){
        for(i=0; i<this->system->getNParticles(); i++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
            if (j==0) this->system->getParticles()[i]->setPosition(pos_var);
        }
        
    }
    */

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }

        this->system->getParticles()[i]->setPosition(pos_var);
    }
    
    if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}

    for(i=1; i<thermalization_step; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));
        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
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
        
    }


    for(i=thermalization_step; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));
        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
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

        if( i>=(int)(thermalization_step) ){
            
            if(i==(int)(thermalization_step)){
                tmp1 = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
            } else {
                if(last_accepted){
                    tmp1 = this->system->getHamiltonian()->LocalEnergyAnalytic();
                }
            }
            energy += tmp1;
            energy2 += tmp1*tmp1;

            if (allAverages){
                tmp2 = this->system->getWavefunction()->psibar_psi();
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


vector<double> Metropolis::solve(double r_max, int N_bins){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, k=0, idx=0;
    double energy=0.0, energy2=0.0, tmp1=0.0;
    double psi_old = 0.0, psi_new=0.0;
    double dist=0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    bool last_accepted;
    int thermalization_step = (int) this->Nsteps*this->InitialFraction;
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);
    vector<int> counts(N_bins, 0);
    vector<double> r(N_bins+1, 0.0);
    for(i=0; i<N_bins; i++){
        r[i] = (double) i* r_max / N_bins;
    }

    /* I DON'T KNOW IF THIS IS WORKING PROPERLY
    #pragma omp parallel for collapse(2)
    for(j=0; j<this->system->getDimension(); j++){
        for(i=0; i<this->system->getNParticles(); i++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
            if (j==0) this->system->getParticles()[i]->setPosition(pos_var);
        }
        
    }
    */

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        this->system->getParticles()[i]->setPosition(pos_var);

        dist = sqrt(this->system->r2(pos_var, (double) 1.0));
        for(k=0; k<(N_bins); k++){
            if((r[k] < dist) && (dist < r[k+1])){
                counts[k]++;
            }
        }
    }
    
    if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}


    for(i=1; i<thermalization_step; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));
        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
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
        
    }


    for(i=thermalization_step; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));
        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
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

        if( i>=(int)(thermalization_step) ){
            
            if(i==(int)(thermalization_step)){
                tmp1 = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
            } else {
                if(last_accepted){
                    tmp1 = this->system->getHamiltonian()->LocalEnergyAnalytic();
                }
            }
            energy += tmp1;
            energy2 += tmp1*tmp1;
            
        }

        for(j=0; j<this->system->getNParticles(); j++){
            dist = sqrt(this->system->r2(this->system->getParticles()[j]->getPosition(), (double) 1.0));
            for(k=0; k<(N_bins); k++){
                if((r[k] < dist) && (dist < r[k+1])){
                    counts[k]++;
                }
            }
        }       
        
    }

    energy = energy/this->Nsteps/(1-this->InitialFraction);
    energy2 = energy2/this->Nsteps/(1-this->InitialFraction); 
    ratio_accepted = (double) accepted/this->Nsteps;

    
    ofstream onebodyFile;
    onebodyFile.open("./plotting/data/onebody_density.csv");
    onebodyFile << "r,counts";
    for(i=0; i<N_bins; i++){
        onebodyFile << endl << r[i] + 0.5 * r_max / N_bins << "," << (double) counts[i] / this->Nsteps/(1-this->InitialFraction) / this->system->getNParticles();
    }
    onebodyFile.close();

    return {energy, energy2 - energy * energy, ratio_accepted};
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
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }
        this->system->getParticles()[i]->setPosition(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}
    }

    //MCsteps
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
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

