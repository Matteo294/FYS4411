#include "functions.h"

Functions::Functions(System* system) { this->system = system;}



vector<vector<double>> Functions::solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas) {
    
    vector<vector<double>> results(Nalphas + 1, vector<double>(3));
    vector<double> results_prov(2, 0.0);
    int i=0, j=0, k=0, idx=0;

    for(k=0; k<=Nalphas; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());

        // set alpha
        this->system->getWavefunction()->setParameter(0, alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas );

        // solve
        results_prov = this->system->getSolver()->solve((bool) 0);
        
        results[k] = {alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas, results_prov[0], results_prov[1]};
        cout << fixed << setprecision(5) << "alpha= " << results[k][0] << "\t energy= " << results[k][1] << "\t std= " << results[k][2] << endl;
    }

    return results;
}


double Functions::gradientDescent(double initialAlpha, double gamma, double tolerance, int Nmax, int Nsteps){

    cout << "Start searching for best alpha..." << endl;
    this->system->getSolver()->setNsteps(Nsteps);
    int i=0;
    double alpha=initialAlpha, deltaAlpha=0.0;
    vector<double> results;

    // first solve --> we generate the first alphanew
    this->system->getWavefunction()->setParameter(0, alpha);
    results = this->system->getSolver()->solve((bool) 1);
    // print current alpha, energy, derivative of energy wrt 
    cout << fixed << setprecision(5) << "iteration=" << i << "\talpha=" << alpha << "\tenergy=" << results[0] << "\tderivative=" << results[2] << "\t|dAlpha|-tol=" << abs(deltaAlpha) - tolerance << endl;
    deltaAlpha = - gamma * results[2];
    alpha = alpha + deltaAlpha;
    i++;

    while((i<Nmax) && (abs(deltaAlpha)>tolerance)){
        this->system->getWavefunction()->setParameter(0, alpha);
        results = this->system->getSolver()->solve((bool) 1);
        cout << fixed << setprecision(8) << "iteration=" << i << "\talpha=" << alpha << "\tenergy=" << results[0] << "\tderivative=" <<  results[2] << "\t|dAlpha|-tol=" << abs(deltaAlpha) - tolerance << endl;
        deltaAlpha = - gamma * results[2];
        alpha = alpha + deltaAlpha;
        i++;
    }

    return alpha;

}