#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include "../System/system.h"
using namespace std;

class Wavefunction{
    
    public:
        Wavefunction(class System* system, int nparams);
        ~Wavefunction();

        // Evaluates the wavefunction in the point in thich the particle are in this moment 
        virtual double evaluateAll() = 0;

        // Evaluates the gaussian contribution relative to the part_idx-th particle.
        virtual double evaluateSing(int part_idx) = 0;
        // Eavluated the second derivative: usefule when the local energy is not analytic
        virtual double numericalSecondDerivative(int part_idx, int direction, double h) = 0; 
        // Evaluates the drift force associated to the part_idx-th particle
        virtual vector<double> DriftForce(int part_idx) = 0;
        // Evaluates derivative of the function with respect to alpha, and divides by the wavefunction
        virtual double psibar_psi() = 0;
        
        // System pointer
        class System* s;

        // Get and set wavefunction parameters
        void setParameter(int idx, double value);
        double getParameter(int idx);
        // Get number of parameters in the chosen wavefunction
        int getNparams();
        
    protected:
        int nparams;
        vector<double> params; // variational parameters
        
};