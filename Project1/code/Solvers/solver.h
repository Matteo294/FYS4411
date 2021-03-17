#pragma once
#include "../System/system.h"
#include <ctime>
#include <vector>
#include <fstream>
using namespace std;

class Solver{
    public:
        ~Solver();
        Solver(class System* system, int Nsteps, int NstepsThermal, int nparams, bool tofile);

    // Getters
        int getNsteps();
        double getNstepsThermal();
        double getParameter(int idx);
        int getnparameter();

    // Setters
        void setNsteps(int Nsteps);
        void setNstepsThermal(int NstepsThermal);
        void setParameter(int idx, double value);
        void setPrintFile(string new_file);

    // Other functions & attributes
        /** Performs a MC simulation to evaluate the energy of the ground state of the system
         *  if allAverages = 1 --> evaluates the derivative of the mean value of the energy with respect to alpha
        **/
        virtual vector<double> solve(bool allAverages) = 0;

        /** Performs a MC simulation to evaluate the energy of the ground state of the system and the onebody density
         * first argument specifies the r_max to sample and second argument is the number of bins to fill
         * the results are printed to file
        **/
        virtual vector<double> solve(double r_max, int N_bins) = 0;

        /** Performs a MC simulation to evaluate the energy of the ground state of the system
         *  This function overrides solve(bool allAverages) fro numerical evaluations. 
         * The parameter h is the steplength used to evaluat numerical derivatives 
         **/
        virtual vector<double> solve(double h) = 0; // override when numerical local energy

        virtual void thermalize() = 0;
        class System* system;

    protected:
        int Nsteps;
        int NstepsThermal;
        int nparams;
        bool tofile;
        vector<double> params;
        ofstream energytofile;
};