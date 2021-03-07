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

        /// Evaluates the wavefunction in the point in thich the particle are in this moment 
        virtual double evaluateAll() = 0;

        /// Evaluates the gaussian contribution relative to the part_idx-th particle. \see evaluateAll()
        virtual double evaluateSing(int part_idx) = 0;

        /** Evaluates numerically the second derivative with respect to the coordinate "direction" of particle "part_idx".
            The derivative is evaluated in the point in which the particles are in this moment.
            direction can be 0 (x), 1 (y), 2 (z), accordingly to the dimension of the system chosen. 
        */
        virtual double numericalSecondDerivative(int part_idx, int direction, double h) = 0; 

        /// Evaluates the drift force associated to the part_idx-th particle
        virtual vector<double> DriftForce(int part_idx) = 0;
        
        /// Evaluates the analytical derivative with respect to alpha
        virtual double analyticalAlphaDerivative() = 0;
        
        class System* s;

        // Get and set variational parameters
        void setParameter(int idx, double value);
        double getParameter(int idx);
        
    protected:
        int nparams;
        vector<double> params; // variational parameters
        
};