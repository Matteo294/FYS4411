#include "functions.h"
#include <omp.h>
#include<string>

Functions::Functions(System* system) { this->system = system;}
Functions::~Functions(){};

vector<double> Functions::solve_singleRun(){

    if(this->system->getSolver()->getToFile()){
        this->system->getSolver()->energytofile = fopen("./Analysis/Data/standard/singlerun/energyateverystep.dat","w");
    }

    this->system->getSolver()->thermalize();
    vector<double> results= this->system->getSolver()->solve(false);
    this->printResultsSolver(results);

    return results;
}

vector<double> Functions::solve_singleRun(double h){
    if(this->system->getSolver()->getToFile()){
        this->system->getSolver()->energytofile = fopen("./Analysis/Data/standard/singlerun/energyateverystep.dat","w");
    }

    this->system->getSolver()->thermalize();
    vector<double> results= this->system->getSolver()->solve(h);
    this->printResultsSolver(results);

    return results;
}


vector<double> Functions::solve_singleRun(double rmax, int Nbins){

    if(this->system->getSolver()->getToFile()){
        this->system->getSolver()->energytofile = fopen("./Analysis/Data/standard/onebody_density/energyateverystep.dat","w");
    }

    this->system->getSolver()->thermalize();
    vector<double> results= this->system->getSolver()->solve(rmax, Nbins);
    this->printResultsSolver(results);

    return results;
}



vector<vector<double>> Functions::solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool alphatoFile) {
    
    if(this->system->getSolver()->getToFile()){
        experimental::filesystem::remove_all("./Analysis/Data/standard/varying_alpha");
        experimental::filesystem::create_directory("./Analysis/Data/standard/varying_alpha");
    }

    
    vector<vector<double>> results(Nalphas + 1);
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;
    double kalpha;

    // print alpha values to csv
    if(alphatoFile){
        this->alphaFile.open("./Analysis/Data/standard/varying_alpha/varying_alpha.csv");
        this->alphaFile << "alpha";
        for(k=0; k<=Nalphas; k++){
            this->alphaFile << endl << alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas;
        }
        this->alphaFile.close();
    }

    for(k=0; k<=Nalphas; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());
        results[k].resize(3);

        kalpha = alpha_min + (double) k * (alpha_max - alpha_min) / Nalphas;
        cout << kalpha << " ";
        //set alpha
        this->system->getWavefunction()->setParameter(0, kalpha);

        if(this->system->getSolver()->getToFile()){
            this->system->getSolver()->setPrintFile("./Analysis/Data/standard/varying_alpha/energyateverystep"+to_string(k));
        }
        
        // solve
        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false);
        results[k]= {kalpha, results_prov[0], results_prov[1], results_prov[2]};
        
        cout << fixed << setprecision(5) << "alpha: " << kalpha << "\t ";
        this->printResultsSolver(results_prov);
        cout << endl;

    }
    

    return results;
}




vector<vector<double>> Functions::solve_varying_dt(double dt_min, double dt_max, int Ndt, bool dttoFile) {
    assert(this->system->getSolver()->getnparameter()==2);
    double expmin = log10(dt_min);
    double expmax = log10(dt_max);
    double expstep = (double) (expmax-expmin)/Ndt;
    vector<vector<double>> results(Ndt + 1, vector<double>(4));
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;
    double kexp, kdt;

    if(this->system->getSolver()->getToFile()){
        experimental::filesystem::remove_all("./Analysis/Data/standard/varying_dt");
        experimental::filesystem::create_directory("./Analysis/Data/standard/varying_dt");
    }

    // print alpha values to csv
    if(dttoFile){
        this->dtFile.open("./Analysis/Data/standard/varying_dt/varying_dt.csv");
        this->dtFile << "dt";
        for(k=0; k<=Ndt; k++){
            this->dtFile << endl << pow(10, expmin + k*expstep);
        }
        this->dtFile.close();
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

        if(this->system->getSolver()->getToFile()){
            this->system->getSolver()->setPrintFile("./Analysis/Data/standard/varying_dt/energyateverystep"+to_string(k));
        }
        
        // solve
        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false);
        results[k] = {kdt, results_prov[0], results_prov[1], results_prov[2]};
        cout << fixed << setprecision(5) << "dt: " << kdt << "\t";
        this->printResultsSolver(results_prov);
        cout << endl;
    }


    return results;
}


vector<vector<double>> Functions::solve_varying_N(vector<int> N, bool NtoFile){
    assert(this->system->getNParticles() < N.at(0) );
    vector<vector<double>> results(N.size());
    vector<double> results_prov(3, 0.0);
    vector<double> zeros(this->system->getDimension(), 0.0);
    int n=0, i=0;

    if(this->system->getSolver()->getToFile()){
        experimental::filesystem::remove_all("./Analysis/Data/standard/varying_N");
        experimental::filesystem::create_directory("./Analysis/Data/standard/varying_N");
    }

    if(NtoFile){
        this->NFile.open("./Analysis/Data/standard/varying_N/varying_N.csv");
        this->NFile << "N";

        for(n=0; n<N.size(); n++){
            this->NFile << endl << N[n];
        }
        this->NFile.close();

    }

    for(n=0; n<N.size(); n++){
        i = this->system->getNParticles();
        while(i<N[n]){
            this->system->addParticle(1.0, zeros);
            i++;
        }

        if(this->system->getSolver()->getToFile()){
            this->system->getSolver()->setPrintFile("./Analysis/Data/standard/varying_N/energyateverystep"+to_string(n));
        }

        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false); 
        results[n] = { (double) N[n], results_prov[0], results_prov[1], results_prov[2]};
        cout << "N: " << N[n] << "\t "; 
        this->printResultsSolver(results_prov);
        cout << endl;     
    }

    return results;
}



double Functions::gradientDescent(double initialAlpha, double gamma, double tolerance, int NiterMax, int Nsteps){
    
    cout << "Start searching for best alpha..." << endl;
    bool tofile_set=this->system->getSolver()->getToFile();
    //this->system->getSolver().sett
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

    cout << scientific << setprecision(5) << "best alpha= " << alpha << endl;
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