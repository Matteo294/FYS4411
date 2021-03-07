#include "solver.h"

Solver::~Solver(){};

Solver::Solver(System* system, int Nsteps, double initialFraction, int nparams){
    this->system = system;
    this->Nsteps=Nsteps;
    this->InitialFraction=initialFraction;
    this->nparams=nparams;
    this->params.resize(nparams, 0.0);
}


// Getters
    int Solver::getNsteps() { return this->Nsteps; }
    double Solver::getInitialFraction() { return this->InitialFraction; }
    int Solver::getnparameter(){ return this->nparams; }
    double Solver::getParameter(int idx){
        assert(idx < this->nparams);
        return this->params[idx];
    }

// Setters
    void Solver::setNsteps(int Nsteps) { this->Nsteps = Nsteps; }
    void Solver::setInitialFraction(double initialFraction) { this->InitialFraction = initialFraction; }
    void Solver::setParameter(int idx, double value){
        assert(idx < this->nparams);
        this->params[idx] = value;
    }