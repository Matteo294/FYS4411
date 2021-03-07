#include "functions.h"

Functions::Functions(System* system) { this->system = system;}
Functions::~Functions(){};


vector<vector<double>> Functions::solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool toFile) {
    
    vector<vector<double>> results(Nalphas + 1, vector<double>(4));
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;

    if(toFile){
        this->alphaFile.open("./plotting/data/varying_alpha.csv");
        this->alphaFile << "alpha,energy,std";
    }

    for(k=0; k<=Nalphas; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());

        // set alpha
        this->system->getWavefunction()->setParameter(0, alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas );
        
        // solve
        results_prov = this->system->getSolver()->solve((bool) 0);
        
        results[k] = {alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas, results_prov[0], results_prov[1], results_prov[2]};
        cout << fixed << setprecision(5) << "alpha= " << results[k][0] << "\t energy= " << results[k][1] << "\t std= " << results[k][2] << "\taccepted=" << results_prov[2] << endl;
        this->alphaFile << endl << results[k][0] << "," << results[k][1] << "," << results[k][2];
    }
    this->alphaFile.close();

    return results;
}


vector<vector<double>> Functions::solve_varying_dt(double dt_min, double dt_max, int Ndt) {
    assert(this->system->getSolver()->getnparameter()==2);
    vector<vector<double>> results(Ndt + 1, vector<double>(4));
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;

    for(k=0; k<=Ndt; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());

        // set alpha
        this->system->getSolver()->setParameter(0, dt_min + (double) k * (dt_max - dt_min) / Ndt );
        
        // solve
        results_prov = this->system->getSolver()->solve((bool) 0);
        
        results[k] = {dt_min + (double) k * (dt_max - dt_min) / Ndt, results_prov[0], results_prov[1], results_prov[2]};
        cout << fixed << setprecision(5) << "dt= " << results[k][0] << "\t energy= " << results[k][1] << "\t std= " << results[k][2] << "\taccepted=" << results_prov[2] << endl;
    }

    return results;
}


double Functions::gradientDescent(double initialAlpha, double gamma, double tolerance, int Nmax, int Nsteps){
    
    cout << "Start searching for best alpha..." << endl;
    this->system->getSolver()->setNsteps(Nsteps);
    int i=0;
    double alpha=initialAlpha, deltaAlpha=1.0;
    vector<double> results(4, 0.0);

    while((i<Nmax) && (abs(deltaAlpha)>tolerance)){
        this->system->getWavefunction()->setParameter(0, alpha);
        results = this->system->getSolver()->solve(true);
        cout << fixed << setprecision(8) << "iteration=" << i << "\talpha=" << alpha << "\tenergy=" << results[0] << "\tstd=" << results[1] << "\taccepted=" << results[2] << "\tderivative=" << results[3] << endl;
        deltaAlpha = - gamma * results[3];
        alpha = alpha + deltaAlpha;
        i++;
    }

    return alpha;

}
vector<vector<double>> Functions::solve_varying_N(vector<int> N){
    vector<vector<double>> res(N.size());
    for(int n=0; n<N.size(); n++){
        cout << n << " " << this->system->getNParticles() << " " << n - this->system->getNParticles() << endl;
        for(int i=0; i<(N[n] - this->system->getNParticles()); i++){
            this->system->addParticle(1.0, {0.0, 0.0, 0.0});
        }           
        res[n] = this->system->getSolver()->solve(false); 
        cout << this->system->getNParticles() << " ";
    }
    return res;
}

void printToFile(double val, ofstream f, string s) {f << val << s;}