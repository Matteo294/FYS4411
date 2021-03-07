#pragma once
#include "../System/system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *system);
        ~Hamiltonian();
        
        class System* system;

        /** Evaluates the local energy of the system analytically.
         * \attention This analytical form is strictly dependent on the chosen hamiltonian and the chosen wavefunction
         * \see LocalEnergyNumeric(double h)
        **/
        virtual double LocalEnergyAnalytic() = 0;
        /// Evaluates the local energy numerically. Way more computationally expensive than the analytic version.
        /// \see LocalEnergyAnalytic()
        virtual double LocalEnergyNumeric(double h) = 0;

};
