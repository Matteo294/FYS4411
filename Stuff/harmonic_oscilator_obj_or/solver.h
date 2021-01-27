#pragma once
#include "system.h"

class Solver{
    public:
        Solver(class System* system, int Nsteps, double step, double initialFraction);
        virtual int getNsteps()=0;
        virtual double getstep()=0;
        virtual double getinitialFraction()=0;
        virtual void solve()=0;
        int Nsteps;
        double step;
        double initialFraction;
        class System* system;

};