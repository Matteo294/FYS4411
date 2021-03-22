#pragma once
#include "../System/system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *systemm, int nparams);
        ~Hamiltonian();
        
        class System* system;
        void setParameter(int idx, double value);

        double getParameter(int idx);
        int getnparameter();

        // Evaluates the local energy of the system analytically.
        virtual double LocalEnergyAnalytic() = 0;
        // Evaluates the local energy numerically
        virtual double LocalEnergyNumeric(double h) = 0;

    protected:
        int nparams;
        vector<double> params;

};
