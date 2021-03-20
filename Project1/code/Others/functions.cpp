#include "functions.h"
#include <omp.h>
#include<string>

Functions::Functions(System* system, bool parallel) { 
    this->system = system;
    this->parallel = parallel;
}
Functions::~Functions(){};

vector<double> Functions::solve_singleRun(){

    if(this->system->getSolver()->getToFile()){

        if(this->parallel){
            this->system->getSolver()->setPrintFile("./Analysis/Data/parallel/singlerun/energyateverystep" + to_string(omp_get_thread_num()));
        } else {
            this->system->getSolver()->energytofile = fopen("./Analysis/Data/standard/singlerun/energyateverystep.dat","w");
        }
        
    }

    //cout << "qui" << endl;
    this->system->getSolver()->thermalize();
    vector<double> results= this->system->getSolver()->solve(false);
    this->printResultsSolver(results);
    

    return results;
}

vector<double> Functions::solve_singleRun(double h){
    if(this->system->getSolver()->getToFile()){
        if(this->parallel){
            this->system->getSolver()->setPrintFile("./Analysis/Data/parallel/singlerun/energyateverystep" + to_string(omp_get_thread_num()));
        } else {
            this->system->getSolver()->energytofile = fopen("./Analysis/Data/standard/singlerun/energyateverystep.dat","w");
        }
    }

    this->system->getSolver()->thermalize();
    vector<double> results= this->system->getSolver()->solve(h);
    this->printResultsSolver(results);

    return results;
}




vector<vector<double>> Functions::solve_varying_alpha(double alpha_min, double alpha_max, int Nalphas, bool alphatoFile) {
    
    vector<vector<double>> results(Nalphas + 1);
    vector<double> results_prov(3, 0.0);
    int i=0, j=0, k=0, idx=0;
    double kalpha;

    // print alpha values to csv
    if(alphatoFile){
        if( (this->parallel && (omp_get_thread_num()==0)) ){
            this->alphaFile.open("./Analysis/Data/parallel/varying_alpha/varying_alpha.csv");
        } else if (!this->parallel ){
            this->alphaFile.open("./Analysis/Data/standard/varying_alpha/varying_alpha.csv");
        }

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
        //set alpha
        this->system->getWavefunction()->setParameter(0, kalpha);

        if(this->system->getSolver()->getToFile()){
            if(this->parallel){
                this->system->getSolver()->setPrintFile("./Analysis/Data/parallel/varying_alpha/energyateverystepcore" + to_string(omp_get_thread_num()) + "alpha" + to_string(k));
            } else {
                this->system->getSolver()->setPrintFile("./Analysis/Data/standard/varying_alpha/energyateverystepalpha" + to_string(k));    
            }
        }
        
        // solve
        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false);
        results[k]= {kalpha, results_prov[0], results_prov[1], results_prov[2]};
        
        if( (this->parallel && (omp_get_thread_num()==0)) || (!this->parallel)){
            cout << fixed << setprecision(5) << "alpha: " << kalpha << "\t ";
        }
        
        this->printResultsSolver(results_prov);

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



    // print dt values to csv
    if(dttoFile){
        if( (this->parallel && (omp_get_thread_num()==0)) ){
            this->dtFile.open("./Analysis/Data/parallel/varying_dt/varying_dt.csv");
            this->dtFile << "dt";
            for(k=0; k<=Ndt; k++){
                this->dtFile << endl <<  pow(10, expmin + k*expstep);

            }
            this->dtFile.close();
        } else if (!this->parallel ){
            this->dtFile.open("./Analysis/Data/standard/varying_dt/varying_dt.csv");
            this->dtFile << "dt,acceptance";
        }

        
    }

    for(k=0; k<=Ndt; k++){
        // initialize random variable
        random_device rd;
        mt19937_64 gen(rd());

        kexp = expmin + k*expstep;
        kdt = pow(10,kexp);
        //cout << kdt << " " << endl;

        // set dt
        this->system->getSolver()->setParameter(0, kdt);

        if(this->system->getSolver()->getToFile()){
            if(this->parallel){
                this->system->getSolver()->setPrintFile("./Analysis/Data/parallel/varying_dt/energyateverystepcore" + to_string(omp_get_thread_num()) + "dt" + to_string(k));
            } else {
                this->system->getSolver()->setPrintFile("./Analysis/Data/standard/varying_dt/energyateverystepdt" + to_string(k));    
            }
        }
        
        // solve
        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false);
        results[k] = {kdt, results_prov[0], results_prov[1], results_prov[2]};

        if(dttoFile && !this->parallel){
            this->dtFile << endl <<  kdt << "," << results_prov[2];
        }

        if( (this->parallel && (omp_get_thread_num()==0)) || (!this->parallel)){
            cout << fixed << setprecision(5) << "dt: " << kdt << "\t ";
        }
        
        this->printResultsSolver(results_prov);
    }
    
    if(dttoFile && !this->parallel)
    {
        this->dtFile.close();
    }


    return results;
}


vector<vector<double>> Functions::solve_varying_N(vector<int> N, bool NtoFile){
    assert(this->system->getNParticles() <= N.at(0) );
    vector<vector<double>> results(N.size());
    vector<double> results_prov(3, 0.0);
    vector<double> zeros(this->system->getDimension(), 0.0);
    int n=0, i=0;

    if(NtoFile){
        if( (this->parallel && (omp_get_thread_num()==0)) ){
            this->NFile.open("./Analysis/Data/parallel/varying_N/varying_N.csv");
        } else if (!this->parallel ){
            this->NFile.open("./Analysis/Data/standard/varying_N/varying_N.csv");
        }

        this->NFile << "N";
        for(n=0; n<=N.size(); n++){
            this->NFile << endl <<  N[n];
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
            if(this->parallel){
                this->system->getSolver()->setPrintFile("./Analysis/Data/parallel/varying_N/energyateverystepcore" + to_string(omp_get_thread_num()) + "N" + to_string(n));
            } else {
                this->system->getSolver()->setPrintFile("./Analysis/Data/standard/varying_N/energyateverystepN" + to_string(n));    
            }
        }


        this->system->getSolver()->thermalize();
        results_prov = this->system->getSolver()->solve(false); 
        results[n] = { (double) N[n], results_prov[0], results_prov[1], results_prov[2]};

        
        if( (this->parallel && (omp_get_thread_num()==0)) || (!this->parallel)){
            cout << fixed << setprecision(5) << "N: " << N[n] << "\t ";
        }
        
        this->printResultsSolver(results_prov);  
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



vector<double> Functions::solve_singleRun(double rmax, int Nbins){

    if(this->system->getSolver()->getToFile()){

        if(this->parallel){
            this->system->getSolver()->setPrintFile("./Analysis/Data/parallel/onebody_density/energyateverystep" + to_string(omp_get_thread_num()));
        } else {
            this->system->getSolver()->setPrintFile("./Analysis/Data/standard/onebody_density/energyateverystep");
        }
        
    }

    if(this->parallel){
         this->system->getSolver()->setOneBodyFile("./Analysis/Data/parallel/onebody_density/counts" + to_string(omp_get_thread_num()));
    } else {
        this->system->getSolver()->setOneBodyFile("./Analysis/Data/standard/onebody_density/counts");
    }
    

    this->system->getSolver()->thermalize();
    vector<double> results = this->system->getSolver()->solve(rmax, Nbins);
    this->printResultsSolver(results);
    

    return results;
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
    if(this->parallel && (omp_get_thread_num()==0)){
        cout << scientific << setprecision(5) << "core#0--> E: " << res[0] << "\t var: " << res[1] << fixed << "\t acceptance: " << res[2] << " " << endl;
    } else if(!this->parallel) {
        cout << scientific << setprecision(5) << "E: " << res[0] << "\t var: " << res[1] << fixed << "\t acceptance: " << res[2] << endl;
    }
}