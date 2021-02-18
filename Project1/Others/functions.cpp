#include "functions.h"

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

void Functions::bestAlpha(double gamma){
    double avg_psidot_psi_eloc = this->system->getSolver()->computeExpectationPsidotPsiEl()[0];
    double avg_psidot_psi = this->system->getSolver()->computeExpectationPsidotPsi()[0];
    double avg_energy = this->system->getSolver()->computeEnergy()[0];

    double alpha_old=1.0, alpha_new=0.51; // just to make sure to enter the cycle
    int i=0;
    while((abs(alpha_new - alpha_old) > 1e-6) && (i<1000000)){
        alpha_old = alpha_new;

        system->getWavefunction()->setParameter(0, alpha_new);

        avg_energy = this->system->getSolver()->computeEnergy()[0];
        avg_psidot_psi = this->system->getSolver()->computeExpectationPsidotPsi()[0];
        avg_psidot_psi_eloc = this->system->getSolver()->computeExpectationPsidotPsiEl()[0];

        alpha_new = alpha_old - 2*gamma*(avg_psidot_psi_eloc - avg_psidot_psi*avg_energy);
        cout << alpha_new << " " << avg_energy << " " << 2*gamma*(avg_psidot_psi_eloc - avg_psidot_psi*avg_energy) << endl;
        i++;
    }
    cout << alpha_new << endl;
}