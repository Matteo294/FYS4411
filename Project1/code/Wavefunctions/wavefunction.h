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

        class System* s;

        /** Evaluates the wavefunction in the point in which the particle are in this moment 
        **/
        virtual double evaluateAll() = 0;
        
        /** Evaluates the terms of the wavefunctions that contain information relative to the part_idx-th particle.
         * **/
        virtual double evaluateSing(int part_idx) = 0;
        
        /** Evaluates numerically the second derivative of the terms of the wavefunction containing information on part_idx-th particle.
        * This is usefull when the analytic expression for the local energy is not known.
        **/
        virtual double numericalSecondDerivative(int part_idx, int direction, double h) = 0; 
        
        /** Evaluates the drift force associated to the part_idx-th particle
        **/
        virtual vector<double> DriftForce(int part_idx) = 0;
        
        /** Evaluates derivative of the function with respect to alpha, and divides by the wavefunction
        **/
        virtual double psibar_psi() = 0;

        /** Set wavefunction parameters in an array
        **/
        void setParameter(int idx, double value);
        
        /** Get wavefunction parameters in an array
        **/
        double getParameter(int idx);
        
        /** Get number of parameters in the chosen wavefunction
        **/
        int getNparams();
        
    protected:
        int nparams; // dimension of the array containing the parameters
        vector<double> params; // variational parameters
        
};