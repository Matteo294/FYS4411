#include "solver.h"

Solver::Solver(System* system, int Nsteps, double initialFraction){
    this->system = system;
    this->Nsteps=Nsteps;
    this->InitialFraction=initialFraction;
}


// Getters
    int Solver::getNsteps() { return this->Nsteps; }
    double Solver::getInitialFraction() { return this->InitialFraction; }

// Setters
    void Solver::setNsteps(int Nsteps) { this->Nsteps = Nsteps; }
    void Solver::setInitialFraction(double initialFraction) { this->InitialFraction = initialFraction; }