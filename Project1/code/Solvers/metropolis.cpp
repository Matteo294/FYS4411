#include "metropolis.h"
#include <omp.h>

Metropolis::~Metropolis(){};

Metropolis::Metropolis(System* system, int Nsteps, int NstepsThermal, double step, bool tofile) : Solver(system, Nsteps, NstepsThermal, 1, tofile){
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
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);


    // particles arrive here already thermalized
    for(i=0; i<this->Nsteps; i++){
        
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

        if(i==0 || last_accepted){
            tmp1 = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
        }

        energy += tmp1;
        energy2 += tmp1*tmp1;

        if (allAverages){
            tmp2 = this->system->getWavefunction()->psibar_psi();
            psi_bar_psi += tmp2;
            psi_bar_psi_EL += tmp2 * tmp1;
        }

        if (this->tofile && (i%1==0)){
            fprintf(energytofile, "%f\n",tmp1);
        }
        /*
        if(i%10000==0){
            cout << fixed << setprecision(2) << "\rprogress " << 100 * (double) i / this->Nsteps << "%" << flush;
        }
        */
       
    
    }

    //cout << "\33[2K\r";
    if(this->tofile){ fclose(energytofile); }
    energy = energy/(this->Nsteps);
    energy2 = energy2/(this->Nsteps); 
    psi_bar_psi = psi_bar_psi/this->Nsteps;
    psi_bar_psi_EL = psi_bar_psi_EL/this->Nsteps;
    ratio_accepted = (double) accepted/this->Nsteps;

    return {energy, energy2 - pow(energy, 2), ratio_accepted, 2 * (psi_bar_psi_EL - psi_bar_psi * energy)};
}




vector<double> Metropolis::solve(double h){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp1=0.0;
    double psi_old=0.0, psi_new=0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    bool last_accepted;
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    for(i=0; i<this->Nsteps; i++){
        
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


        if(i==0 || last_accepted){
            tmp1 = (double) this->system->getHamiltonian()->LocalEnergyNumeric(h);
        } 
        energy += tmp1;
        energy2 += tmp1*tmp1;

        if (this->tofile && (i%1==0)){
            fprintf(energytofile, "%f\n",tmp1);
        }

        /*
        if(i%10000==0){
            cout << fixed << setprecision(2) << "\rprogress " << 100 * (double) i / this->Nsteps << "%" << flush;
        }
        */
       
    
    }

    //cout << "\33[2K\r";
    if(this->tofile){ fclose(energytofile); }
    energy = energy/this->Nsteps;
    energy2 = energy2/this->Nsteps; 
    ratio_accepted = (double) accepted/this->Nsteps;
    
    return {energy, energy2 - pow(energy, 2), ratio_accepted};
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
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);
    vector<int> counts(N_bins, 0);
    vector<double> r(N_bins+1, 0.0);
    for(i=0; i<N_bins; i++){
        r[i] = (double) i* r_max / N_bins;
    }

    for(i=0; i<this->Nsteps; i++){
        
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

        
        if(i==0 || last_accepted){
            tmp1 = (double) this->system->getHamiltonian()->LocalEnergyAnalytic();
        }

        for(j=0; j<this->system->getNParticles(); j++){
            dist = sqrt(this->system->r2(this->system->getParticles()[j]->getPosition(), (double) 1.0));
            for(k=0; k<(N_bins); k++){
                if((r[k] < dist) && (dist < r[k+1])){
                    counts[k]++;
                }
            }
        }
        
        energy += tmp1;
        energy2 += tmp1*tmp1;

        if (this->tofile && (i%1==0)){
            fprintf(energytofile, "%f\n",tmp1);
        }

        /*
        if(i%10000==0){
            cout << fixed << setprecision(2) << "\rprogress " << 100 * (double) i / this->Nsteps << "%" << flush;
        }
        */

    }

    //cout << "\33[2K\r";
    if(this->tofile){ fclose(energytofile); }
    energy = energy/this->Nsteps;
    energy2 = energy2/this->Nsteps; 
    ratio_accepted = (double) accepted/this->Nsteps;

    
    fprintf(this->onebodyFile,"r,counts\n");

    for(i=0; i<N_bins; i++){
        fprintf(this->onebodyFile, "%f,%f\n", r[i] + 0.5 * r_max / N_bins, (double) counts[i] / this->Nsteps/1 / this->system->getNParticles());
    }

    fclose(this->onebodyFile);

    return {energy, energy2 - energy * energy, ratio_accepted};

}






void Metropolis::thermalize(){
    // initialize random variable
    random_device rd;
    mt19937_64 gen(rd());
    
    int i=0, j=0, idx=0;
    double psi_old = 0.0, psi_new=0.0;
    bool usematrix=this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = 2*this->params[0]*( this->system->getRandomGenerator()->uniform(gen) - 0.5);
        }

        this->system->getParticles()[i]->setPosition(pos_var);
    }
    
    if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}

    for(i=0; i<this->NstepsThermal; i++){
        
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
        } 

        /*
        if(i%1000==0){
            cout << fixed << setprecision(2) << "\rprogress thermalization " << 100 * (double) i / this->NstepsThermal << "%" << flush;
        }
        */

    }
    //cout << "\33[2K\r";
    
}