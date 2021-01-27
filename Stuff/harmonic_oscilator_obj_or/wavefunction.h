#pragma once
#include <vector>
#include "system.h"

/* This is the wavefunction superclass. It contains the general features of the wf
such as a function that evaluates the wf in a particular point or a fuction that evaluate
the second derivative in a particular point. */

class Wavefunction{
    
    public:
        Wavefunction(class System* system);

        /* evaluate() evaluates the wf in the specified point. 
        "Virtual" means that it must be implemented in the subclass, 
        for the moment we just create it let others know that it must be implemented */
        virtual double evaluate() = 0;

        virtual double evaluateSecondDerivative() = 0; // Evaluates the second derivative if an analytical form exists
        virtual double numericalSecondDerivative() = 0; // Evaluates the second derivative numerically from the analytical wf
        
        /* This is a pointer to the system. We could have passed the system entirely but it would have slown down the code a lot.
        Just think that for 10 particles and 10000 integration steps we end up with 6x10x10000 variables (positions + velocities) that 
        should be passed one by one if we pass the entire system. Instead by passing the pointer we pass one single variable, the memory address */
        class System* s;

};