#pragma once
#include "../System/system.h"

class Hamiltonian{
    public:
        Hamiltonian(class System *systemm, int nparams);
        ~Hamiltonian();
        
        class System* system;
        /** Assigns the value of the idx-th element of the vector of parameters. 
         * \note The meaning of each element of the vector depends on the chosen Hamiltonian.
         * \see Gaussian
         * \see AsymmeGaussian **/
        void setParameter(int idx, double value);

        /** Returns the idx-th parameter characterizing the chosen Hamiltonian.
         * \note The meaning of each element of the vector depends on the chosen Hamiltonian.
         * \see Gaussian
         * \see AsymmeGaussian **/ 
        double getParameter(int idx);

        /** Returns the number of parameters characterizing the selected Hamiltonian. **/
        int getnparameter();

        /** Evaluates the local energy of the system analytically. **/
        virtual double LocalEnergyAnalytic() = 0;
        /** Evaluates the local energy of the system analytically. **/
        virtual double LocalEnergyNumeric(double h) = 0;

    protected:
        int nparams;
        vector<double> params;

};
