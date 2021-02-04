#include "solver.h"

Solver::Solver(System* system, int Nsteps, double step, double initialFraction){
    this->system = system;
    this->Nsteps=Nsteps;
    this->step=step;
    this->initialFraction=initialFraction;
}