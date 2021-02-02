#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include "system.h"
using namespace std;

/* This is the wavefunction superclass. It contains the general features of the wf
such as a function that evaluates the wf in a particular point or a fuction that evaluate
the second derivative in a particular point. */

class Wavefunction{
    
    public:
        Wavefunction(class System* system, int nparams);

        /* evaluate() evaluates the wf in the specified point. 
        "Virtual" means that it must be implemented in the subclass, 
        for the moment we just create it let others know that it must be implemented */
        virtual double evaluateAll() = 0;

        /* evaluateSing( int part_idx) evaluates the gaussian contribution relative to the part_idx-th particle */
        virtual double evaluateSing(int part_idx) = 0;

        virtual double evaluateSecondDerivative() = 0; // Evaluates the second derivative if an analytical form exists

        /* Evaluates numerically the second derivative with respect to the coordinate "direction" of particle "part_idx".
         The derivative is evaluated in the point in which the particles are in this moment.
        direction can be 0 (x), 1 (y), 2 (z), accordingly to the dimension of the system chosen. 
        */
        virtual double numericalSecondDerivative(int part_idx, int direction, double h) = 0; 

        /* evaluates the drift force associated to the part_idx-th particle */
        virtual vector<double> DriftForce(int part_idx) = 0;
        
        /* This is a pointer to the system. We could have passed the system entirely but it would have slown down the code a lot.
        Just think that for 10 particles and 10000 integration steps we end up with 6x10x10000 variables (positions + velocities) that 
        should be passed one by one if we pass the entire system. Instead by passing the pointer we pass one single variable, the memory address */
        class System* s;

        // Get and set variational parameters
        double setParameter(int idx, double value);
        double getParameter(int idx);
        int nparams;
        std::vector<double> params; // variational parameters
        
};