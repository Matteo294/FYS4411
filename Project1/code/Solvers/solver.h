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
        bool getToFile();


        // Setters
        void setNsteps(int Nsteps);
        void setNstepsThermal(int NstepsThermal);
        void setParameter(int idx, double value);
        void setPrintFile(string new_file);
        void setOneBodyFile(string new_file);

        // Other functions & attributes
        /* Performs a MC simulation to evaluate the energy of the ground state of the system
         *  if allAverages = 1 --> evaluates the derivative of the mean value of the energy with respect to alpha
         */
        virtual vector<double> solve(bool allAverages) = 0;

        /* Performs a MC simulation to evaluate the energy of the ground state of the system and the onebody density
         * first argument specifies the r_max to sample and second argument is the number of bins to fill
         * the results are printed to file
         */
        virtual vector<double> solve(double r_max, int N_bins) = 0;

        /* Performs a MC simulation to evaluate the energy of the ground state of the system
         *  This function overrides solve(bool allAverages) fro numerical evaluations. 
         * The parameter h is the steplength used to evaluat numerical derivatives 
         */
        virtual vector<double> solve(double h) = 0; // override when numerical local energy

        // Thermalizes the system. The number of thermalization steps can be set through the proper setter
        virtual void thermalize() = 0;

        // System pointer
        class System* system;

        // Files to write data
        FILE* energytofile;
        FILE* onebodyFile;

    protected:
        int Nsteps;
        int NstepsThermal;
        int nparams;
        bool tofile;
        vector<double> params;
};