#include "importanceSampling.h"
#include <fstream>

ImportanceSampling::~ImportanceSampling(){};


ImportanceSampling::ImportanceSampling(System* system, int Nsteps, double initialFraction, double dt, double D, bool tofile=0) : Solver(system, Nsteps, initialFraction, 2, tofile){
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

    // particles arrive here already thermalized
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->params[1] * drift_old[j] * this->params[0] + this->system->getRandomGenerator()->normal(gen) * sqrt(this->params[0]);
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

        if (this->tofile){
            fprintf(energytofile, "%f\n",tmp1);
        }
        
        if(i%1000==0 && ( (this->system->getParallel() && omp_get_thread_num()==0) || !this->system->getParallel() )){
            cout << fixed << setprecision(2) << "\rprogress " << 100 * (double) i / this->Nsteps << "%" << flush;
        }
        
    
    }
    
    if(tofile) fclose(this->energytofile);
    cout << "\33[2K\r";
    energy = energy/this->Nsteps;
    energy2 = energy2/this->Nsteps;
    psi_bar_psi = psi_bar_psi/this->Nsteps;
    psi_bar_psi_EL = psi_bar_psi_EL/this->Nsteps;
    ratio_accepted = (double) accepted/this->Nsteps;

    return {energy, energy2 - pow(energy, 2), ratio_accepted, 2 * (psi_bar_psi_EL - psi_bar_psi * energy)};
}





vector<double> ImportanceSampling::solve(double h){
    random_device rd;
    mt19937_64 gen(rd());
        
    int i=0, j=0, idx=0;
    double energy=0.0, energy2=0.0, tmp1=0.0;
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

    // particles arrive here already thermalized
    for(i=1; i<=this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->params[1] * drift_old[j] * this->params[0] + this->system->getRandomGenerator()->normal(gen) * sqrt(this->params[0]);
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

        if(i==0 || last_accepted){
            tmp1 = (double) this->system->getHamiltonian()->LocalEnergyNumeric(h);
        } 

        energy += tmp1;
        energy2 += tmp1*tmp1;

        if (this->tofile){
            fprintf(energytofile, "%f\n", tmp1);
        }

        if(i%1000==0 && ( (this->system->getParallel() && omp_get_thread_num()==0) || !this->system->getParallel() )){
            cout << fixed << setprecision(2) << "\rprogress " << 100 * (double) i / this->Nsteps << "%" << flush;
        }
    }

    if(tofile) fclose(this->energytofile);
    cout << "\33[2K\r";
    energy = energy/this->Nsteps;
    energy2 = energy2/this->Nsteps;
    ratio_accepted = (double) accepted/this->Nsteps;
    
    return {energy, energy2 - pow(energy, 2), ratio_accepted};
}



vector<double> ImportanceSampling::solve(double r_max, int N_bins){
    random_device rd;
    mt19937_64 gen(rd());
        
    int i=0, j=0, k=0, idx=0;
    double energy=0.0, energy2=0.0, tmp1=0.0;
    double psi_old = 0.0, psi_new=0.0;
    double arg = 0.0;
    double dist=0.0;
    int accepted=0;
    double ratio_accepted=0.0;
    bool last_accepted;
    bool usematrix = this->system->getUseMatrix();
    vector<double> pos_old(this->system->getDimension(), 0.0);
    vector<double> pos_var(this->system->getDimension(), 0.0);

    vector<double> drift_old(this->system->getDimension(), 0.0);
    vector<double> drift_new(this->system->getDimension(), 0.0);

    vector<int> counts(N_bins, 0);
    vector<double> r(N_bins+1, 0.0);
    for(i=0; i<N_bins; i++){
        r[i] = (double) i * r_max / N_bins;
    }

    // particles arrive here already thermalized
    for(i=0; i<this->Nsteps; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->params[1] * drift_old[j] * this->params[0] + this->system->getRandomGenerator()->normal(gen) * sqrt(this->params[0]);
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

        if(i==0 || last_accepted){
            tmp1 = this->system->getHamiltonian()->LocalEnergyAnalytic();
        }

        energy += tmp1;
        energy2 += tmp1*tmp1;

        for(j=0; j<this->system->getNParticles(); j++){
            dist = sqrt(this->system->r2(this->system->getParticles()[j]->getPosition(), (double) 1.0));
            for(k=0; k<(N_bins); k++){
                if((r[k] < dist) && (dist < r[k+1])){
                    counts[k]++;
                }
            }
        }

        if (this->tofile){
            fprintf(energytofile, "%f\n", tmp1);
        }

        if(i%1000==0 && ( (this->system->getParallel() && omp_get_thread_num()==0) || !this->system->getParallel() )){
            cout << fixed << setprecision(2) << "\rprogress " << 100 * (double) i / this->Nsteps << "%" << flush;
        }
    
    }

    cout << "\33[2K\r";
    if(this->tofile){ fclose(energytofile); }
    energy = energy/this->Nsteps;
    energy2 = energy2/this->Nsteps;
    ratio_accepted = (double) accepted/this->Nsteps;

    fprintf(this->onebodyFile,"r,counts\n");
    for(i=0; i<N_bins; i++){
        fprintf(this->onebodyFile, "%f,%f\n", r[i], (double) counts[i] / this->Nsteps / this->system->getNParticles());
    }
    
    fclose(this->onebodyFile);
    return {energy, energy2 - pow(energy,2), ratio_accepted};
}


void ImportanceSampling::thermalize(){
    random_device rd;
    mt19937_64 gen(rd());
        
    int i=0, j=0, idx=0;
    double psi_old = 0.0, psi_new=0.0;
    double arg = 0.0;
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
    }

    if(usematrix){ this->system->EvaluateRelativePosition(); this->system->EvaluateRelativeDistance();}

    for(i=1; i<=this->NstepsThermal; i++){
        
        idx = (int) round( this->system->getRandomGenerator()->uniform(gen) * (this->system->getNParticles() - 1));

        pos_old = this->system->getParticles()[idx]->getPosition();
        psi_old = this->system->getWavefunction()->evaluateSing(idx);
        drift_old = this->system->getWavefunction()->DriftForce(idx);

        for(j=0; j<this->system->getDimension(); j++){
            pos_var[j] = this->params[1] * drift_old[j] * this->params[0] + this->system->getRandomGenerator()->normal(gen) * sqrt(this->params[0]);
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
        } 

        if(i%1000==0 && ( (this->system->getParallel() && omp_get_thread_num()==0) || !this->system->getParallel() )){
            cout << fixed << setprecision(2) << "\rthermalization progress " << (double) 100 * i / this->NstepsThermal << "%" << flush;
        }       
       
    }

    cout << "\33[2K\r";
}