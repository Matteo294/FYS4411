#pragma once
#include "solver.h"
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;



class ImportanceSampling : public Solver{
    public:
        ImportanceSampling(class System* system, int Nsteps, double step, double initialFraction, double dt, double D);
        int getNsteps();
        double getstep();
        double getinitialFraction();
        double getdt();
        double getD();
        vector<double> solve(); 
    
    private:
        double dt;
        double D;


};