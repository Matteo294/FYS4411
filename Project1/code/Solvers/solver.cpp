#include "solver.h"
#include <string>
#include <string.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>

Solver::~Solver(){};

Solver::Solver(System* system, int Nsteps, int NstepsThermal, int nparams, bool tofile){
    this->system = system;
    this->Nsteps=Nsteps;
    this->NstepsThermal=NstepsThermal;
    this->nparams=nparams;
    this->tofile = tofile;
    this->params.resize(nparams, 0.0);
    if (tofile){
        energytofile=fopen("./variance/energyateverystep.dat","w");
    }
}

// Getters
    int Solver::getNsteps() { return this->Nsteps; }
    double Solver::getNstepsThermal() { return this->NstepsThermal; }
    int Solver::getnparameter(){ return this->nparams; }
    double Solver::getParameter(int idx){
        assert(idx < this->nparams);
        return this->params[idx];
    }

// Setters
    void Solver::setNsteps(int Nsteps) { this->Nsteps = Nsteps; }
    void Solver::setNstepsThermal(int NstepsThermal) { this->NstepsThermal = NstepsThermal; }
    void Solver::setParameter(int idx, double value){
        assert(idx < this->nparams);
        this->params[idx] = value;
    }
    
    void Solver::setPrintFile(string new_file){
        fclose(energytofile);
        char* char_arr;
        string str_obj("./variance"+ new_file+".dat");
        char_arr = &str_obj[0];
        fopen(char_arr,"w");
    }