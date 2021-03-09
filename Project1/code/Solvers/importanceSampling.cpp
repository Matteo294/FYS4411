#include "importanceSampling.h"
#include <fstream>

ImportanceSampling::~ImportanceSampling(){};


ImportanceSampling::ImportanceSampling(System* system, int Nsteps, double initialFraction, double dt, double D) : Solver(system, Nsteps, initialFraction, 2){
    this->setParameter(0, dt);
    this->setParameter(1, D);
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
    bool usematrix = this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    vector<double> drift_old(this->system->getDimension(), 0.0);
    vector<double> drift_new(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->params[0]);
        }

        this->system->getParticles()[i]->setPosition(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}
    }

    //MCsteps
    for(i=1; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->params[1] * drift_old[j] * this->params[0] + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getParameter(0));
        }

        this->system->getParticles()[idx]->move(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
        psi_new = this->system->getWavefunction()->evaluateSing(idx);
        drift_new = this->system->getWavefunction()->DriftForce(idx);
        
        arg = 0.0;
        for(j=0; j<this->system->getDimension(); j++){
            arg += (drift_old[j] + drift_new[j]) * ( -2 * pos_var[j] + this->params[1] * this->params[0] * (drift_old[j] - drift_new[j]) );
        }

        arg *= 0.25;

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(psi_new, 2) / pow(psi_old, 2))){
            this->system->getParticles()[idx]->setPosition(pos_old);
            if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
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


/* here there are a LOT of things to fix, but we need to ask to morten before proceeding*/
vector<double> ImportanceSampling::solve(double h, bool tofile=0){
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp=0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    double psi_old = 0.0, psi_new = 0.0;
    double arg = 0.0;
    bool usematrix= this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);
    vector<double> drift_old(this->system->getDimension(), 0.0);
    vector<double> drift_new(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->system->getRandomGenerator()->normal(gen) * sqrt(this->params[0]);
        }

        this->system->getParticles()[i]->setPosition(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}
    }



    //MCsteps
    for(i=1; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->params[1] * drift_old[j] * this->params[0] + this->system->getRandomGenerator()->normal(gen) * sqrt(this->getParameter(0));
        }
        
        this->system->getParticles()[idx]->move(pos_var);
        if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
        psi_new = this->system->getWavefunction()->evaluateSing(idx);
        drift_new = this->system->getWavefunction()->DriftForce(idx);
        
        for(j=0; j<this->system->getDimension(); j++){
            arg += (drift_old[j] + drift_new[j]) * ( -2 * pos_var[j] + this->params[1] * this->params[0] * (drift_old[j] - drift_new[j]) );
        }

        arg *= 0.25;
        

        if( this->system->getRandomGenerator()->uniform(gen) > (exp(arg) * pow(psi_new, 2) / pow(psi_old, 2))){
            this->system->getParticles()[idx]->setPosition(pos_old);
            if(usematrix){ this->system->EvaluateRelativePosition(idx); this->system->EvaluateRelativeDistance(idx);}
        } else {
            accepted++;
        }

        if(i>=(int)(this->Nsteps*this->InitialFraction)){
            tmp = (double) this->system->getHamiltonian()->LocalEnergyNumeric(h);
            energy += tmp;
            energy2 += tmp*tmp;
            if(tofile){
                energytofile << tmp << endl; 
            }
        }
        
    }
    energytofile.close();
    energy = energy/this->Nsteps/(1-this->InitialFraction);
    energy2 = energy2/this->Nsteps/(1-this->InitialFraction);
    ratio_accepted = (double) accepted/this->Nsteps;

    return {energy, energy2 - pow(energy, 2), ratio_accepted};
}


