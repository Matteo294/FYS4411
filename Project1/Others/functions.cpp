#include "functions.h"
#include <fstream>

Functions::Functions(System* system) { this->system = system;}




vector< vector<double> > Functions::solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas) {
    
    vector< vector<double> > results(Nalphas + 1, vector<double>(3));
    vector<double> results_prov(2, 0.0);
    int i=0, j=0, k=0, idx=0;

    for(k=0; k<=Nalphas; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());

        // set alpha
        this->system->getWavefunction()->setParameter(0, alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas );

        // solve
        results_prov = this->system->getSolver()->computeEnergy();
        
        results[k] = {alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas, results_prov[0], results_prov[1]};
        cout << fixed << setprecision(5) << "alpha= " << results[k][0] << "\t energy= " << results[k][1] << "\t std= " << results[k][2] << endl;
    }

    return results;
}

void Functions::bestAlpha(double gamma, int nsteps){
    ofstream outputfile; 
    outputfile.open("resultfile.csv");
    outputfile << "alpha" << "," << "averEnergy" << "," << "Eldot" << endl;
    this->system->getSolver()->setNsteps(nsteps);
    double avg_psidot_psi_eloc;
    double avg_psidot_psi;
    double avg_energy;
    double ELdot, old_ELdot;
    double alpha_old=1.0, alpha_new=0.1; // just to make sure to enter the cycle
    int i=0;
    while((abs(alpha_new - alpha_old) > 1e-4) && (i<1000000)){
        alpha_old = alpha_new;

        system->getWavefunction()->setParameter(0, alpha_new);

        avg_energy = this->system->getSolver()->computeEnergy()[0];
        avg_psidot_psi = this->system->getSolver()->computeExpectationPsidotPsi()[0];
        avg_psidot_psi_eloc = this->system->getSolver()->computeExpectationPsidotPsiEl()[0];

        old_ELdot = ELdot;
        ELdot = 2*gamma*(avg_psidot_psi_eloc - avg_psidot_psi*avg_energy);
        alpha_new = alpha_old - ELdot;
        cout << this->system->getWavefunction()->getParameter(0) << " av energy " << avg_energy << "   Eldot  " << ELdot << endl;
        outputfile << this->system->getWavefunction()->getParameter(0) << "," << avg_energy << "," << ELdot << endl;
        i++;
    }
    system->getWavefunction()->setParameter(0, alpha_new);
    cout << "Best alpha " << alpha_new << " found at step " << i << endl;
}