#pragma once
#include "hamiltonian.h"
#include <cmath>

class Elliptical: public Hamiltonian{
    public:

        // Constructor and destructor
        Elliptical(class System* system, double omegaXY, double omegaZ);
        ~Elliptical();

        // Getters
        double getOmegaXY() {return this->omegaXY;}
        double getOmegaZ() {return this->omegaZ;}

        //Other functions
        /** Evaluates the local energy of the system analytically.
         * \attention This analytical form is strictly dependet on the chosen hamiltonian and the chosen wavefunction
         * \see LocalEnergyNumeric(double h), Hamiltonian::LocalEnergyAnalytic()
         **/
        double LocalEnergyAnalytic();
        
        /// \see LocalEnergyAnalytic(), Hamiltonian::LocalEnergyNumeric()
        double LocalEnergyNumeric(double h);

    private:
        // Private attributes: access only via getters and setters
        double omegaXY;
        double omegaZ;
};
