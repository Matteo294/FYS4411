#include "functions.h"
#include <omp.h>

Functions::Functions(System* system) { this->system = system;}
Functions::~Functions(){};


vector<vector<double>> Functions::solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool toFile) {
    
    vector<vector<double>> results(Nalphas + 1);
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;
    double kalpha;

    if(toFile){
        this->alphaFile.open("./plotting/data/varying_alpha.csv");
        this->alphaFile << "alpha,energy,std,acceptance";
    }

    for(k=0; k<=Nalphas; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());
        results[k].resize(3);

        kalpha = alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas;

        // set alpha
        this->system->getWavefunction()->setParameter(0, kalpha);
        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false);
        results[k]= {kalpha, results_prov[0], results_prov[1], results_prov[2]};
        
        // solve
        cout << fixed << setprecision(5) << "alpha: " << kalpha << "\t ";
        this->printResultsSolver(results_prov);
        cout << endl;
        if(toFile) this->alphaFile << endl << results[k][0] << "," << results[k][1] << "," << results[k][2] << "," << results[k][3];
    }
    if(toFile) this->alphaFile.close();

    return results;
}

vector<vector<double>> Functions::solve_varying_dt(double dt_min, double dt_max, int Ndt, bool toFile) {
    assert(this->system->getSolver()->getnparameter()==2);
    double expmin = log10(dt_min);
    double expmax = log10(dt_max);
    double expstep = (double) (expmax-expmin)/Ndt;
    vector<vector<double>> results(Ndt + 1, vector<double>(4));
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;
    double kexp, kdt;

    if(toFile){
        this->dtFile.open("./plotting/data/varying_dt.csv");
        this->dtFile << "dt,energy,std,acceptance";
    }

    for(k=0; k<=Ndt; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());

        kexp = expmin + k*expstep;
        kdt = pow(10,kexp);
        //cout << kdt << " " << endl;

        // set alpha
        this->system->getSolver()->setParameter(0, kdt);
        
        // solve
        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false);
        
        results[k] = {kdt, results_prov[0], results_prov[1], results_prov[2]};
        cout << fixed << setprecision(5) << "dt: " << kdt << "\t";
        this->printResultsSolver(results_prov);
        cout << endl;
        if(toFile) this->dtFile << endl <<  results[k][0] << "," << results[k][1] << "," << results[k][2] << "," << results[k][3];
    }

    if(toFile) this->dtFile.close();

    return results;
}


vector<vector<double>> Functions::solve_varying_N(vector<int> N, bool toFile){
    assert(this->system->getNParticles() < N.at(0) );
    vector<vector<double>> results(N.size());
    vector<double> results_prov(3, 0.0);
    vector<double> zeros(this->system->getDimension(), 0.0);
    int i=0;

    if(toFile){
        this->Nfile.open("./plotting/data/varying_N.csv");
        this->Nfile << "N,energy,std,acceptance";
    }

    for(int n=0; n<N.size(); n++){
        i = this->system->getNParticles();
        while(i<N[n]){
            this->system->addParticle(1.0, zeros);
            i++;
        }
        

        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false); 
        results[n] = { (double) N[n], results_prov[0], results_prov[1], results_prov[2]};
        cout << "N: " << N[n] << "\t "; 
        this->printResultsSolver(results_prov);
        cout << endl;     
        if(toFile) this->Nfile << endl << results[n][0] << "," << results[n][1] << "," << results[n][2] << "," << results[n][3];
    }
    if(toFile) this->Nfile.close();
    return results;
}

double Functions::gradientDescent(double initialAlpha, double gamma, double tolerance, int NiterMax, int Nsteps){
    
    cout << "Start searching for best alpha..." << endl;
    this->system->getSolver()->setNsteps(Nsteps);
    int i=0;
    double alpha=initialAlpha, deltaAlpha=1.0;
    vector<double> results(4, 0.0);

    while((i<NiterMax) && (abs(deltaAlpha)>tolerance)){
        this->system->getWavefunction()->setParameter(0, alpha);
        this->system->getSolver()->thermalize();
        results = this->system->getSolver()->solve(true);
        cout << scientific << setprecision(4) << "iter=" << i << "\talpha=" << alpha << "\tenergy=" << results[0] << "\tstd=" << results[1] << "\tacc=" << results[2] << "\tder_alpha=" << results[3] << endl;
        deltaAlpha = - gamma * results[3];
        alpha = alpha + deltaAlpha;
        i++;
    }
    alpha -= deltaAlpha;
    return alpha;
}





void Functions::printPresentation(){
    cout << endl;
    cout << "***********************************************************" << endl;
    cout << "Emiliano Staffoli, Matteo Zortea, Alexander Ferraro" << endl;
    cout << "Variational Monte Carlo for trapped bosons" << endl;
    cout << "March 2021, University of Oslo" << endl;
    cout << "***********************************************************" << endl;
    cout << endl;
}

void Functions::printResultsSolver(vector<double> res){
    cout << scientific << setprecision(5) << "E: " << res[0] << "\t std: " << res[1] << fixed << "\t acceptance: " << res[2];
}

void Functions::solveParallel(System* s, int N){
    int Nthreads = omp_get_max_threads();
    int Ni = (int) N/Nthreads;
    omp_set_num_threads(Nthreads);
    cout << "Nthreads: " << Nthreads << endl;
    /*omp_set_num_threads(2);
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                s1->getSolver()->setNsteps(Ni);
                cout << s1->getSolver()->solve(false)[0] << endl;
            }
            #pragma omp section
            {
                s2->getSolver()->setNsteps(Ni);
                cout << s2->getSolver()->solve(false)[0] << endl;
            }
        }
    }*/
    #pragma omp parallel for
    for(int i=0; i<Nthreads; i++){
        s->getSolver()->setNsteps(Ni);
        cout << omp_get_thread_num() << "\t" << s->getSolver()->solve(false)[0] << endl;
    } 
}