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
        /** Returns the number of MC steps **/
        int getNsteps();
        /** Returns the number of steps used for thermalizing the system **/
        double getNstepsThermal();
        /** Returns the needed parameters for Importance sampling or Metropolis.
         * \note Their actual meaning depends on the specific choice for the solver.**/
        double getParameter(int idx);
        /** Returns the number of parameters **/
        int getnparameter();
        /** Returns True if printing the local energy at every step to file is enabled.**/
        bool getToFile();


        // Setters
        /** Sets the number of MC steps **/
        void setNsteps(int Nsteps);
        /** Sets the number of steps for thermalizing the system **/
        void setNstepsThermal(int NstepsThermal);
        /** Sets the needed parameters for Importance sampling or Metropolis**/
        void setParameter(int idx, double value);
        /** Sets the file path where to print the local energy**/
        void setPrintFile(string new_file);
        /** Sets the file path where print the local energy **/
        void setOneBodyFile(string new_file);
        /** Sets whether to print to file the local energy at every step. **/
        void setToFile(bool tofile);

        // Other functions & attributes
        /** Performs a MC simulation to evaluate the energy of the ground state of the system
         * if allAverages = 1 --> evaluates the derivative of the mean value of the energy with respect to alpha (used in Functions::gradient_descent)
         **/
        virtual vector<double> solve(bool allAverages) = 0;

        /** Performs a MC simulation to evaluate the energy of the ground state of the system
         *  This function overrides solve(bool allAverages) fro numerical evaluations. 
         * The parameter h is the steplength used to evaluat numerical derivatives 
         **/
        virtual vector<double> solve(double h) = 0; // override when numerical local energy

        /** Performs a MC simulation to evaluate the energy of the system and the radial onebody density
         * first argument specifies the r_max to sample and second argument is the number of bins to fill
         * \note The .csv file with the counts for each bin is always generated. The path can be specified 
         * into Functions::solve_singleRun(double r_max, int N_bins)
         **/
        virtual vector<double> solve(double r_max, int N_bins) = 0;

        /** Thermalizes the system. The number of thermalization steps can be set through the proper setter */
        virtual void thermalize() = 0;

        /** System pointer **/
        class System* system;

        FILE* energytofile;
        FILE* onebodyFile;

    protected:
        int Nsteps;
        int NstepsThermal;
        int nparams;
        bool tofile;
        vector<double> params;
};