#include "config.h"

using namespace std;

int main(int argc, char *argv[]){

const double gamma = 2e-3;

    // Select working mode : run "make" command. Then "./main x" where x is the number indicating working mode (see switch below)
    int selector = 0;
    if(argc>1){
        int mode = stoi(argv[1]);
        assert(mode>=0 && mode<=7);
        selector = stoi(argv[1]);
        if(argc>2){
            activate_asymmetric = (bool) stoi(argv[2]);
            if(argc>3){
                activate_elliptical = (bool) stoi(argv[3]);
                if(argc>4){
                        activate_importance = (bool) stoi(argv[4]);
                }
            }
        }
    }
    
    System system(dimension, Nparticles);

    // Hamiltonians
    Spherical spherical(&system, omegaXY);
    Elliptical elliptical(&system, omegaXY, omegaZ);

    // Wavefunctions
    Gaussian gaussian(&system, alpha);
    AsymmetricGaussian asymmgaussian(&system, alpha, beta, a);

    // Solvers
    Metropolis metropolis(&system, Nsteps_final, NstepsThermal, step, tofile);
    ImportanceSampling importance(&system, Nsteps_final, NstepsThermal, dt, D, tofile);
    
    // Others
    RandomGenerator randomgenerator;
    Functions functions(&system);

    // Choose options (see config.h for changing options and parameters)
    if (activate_elliptical) system.setHamiltonian(&elliptical);
    else system.setHamiltonian(&spherical);

    if (activate_asymmetric) system.setWavefunction(&asymmgaussian);
    else system.setWavefunction(&gaussian);

    if (activate_importance) system.setSolver(&importance);
    else system.setSolver(&metropolis);

    system.setRandomGenerator(&randomgenerator);
    functions.printPresentation();
    
    // !!!!!!!! This should be more precise and should be put inside Function (time should be returned as another value)
    auto start = chrono::steady_clock::now(); // Store starting time to measure run time

    switch(selector){
        case 0: system.getSolver()->thermalize(); 
                functions.printResultsSolver(system.getSolver()->solve(false));  break; // Simple simulation
        case 1: system.getSolver()->thermalize();
                functions.printResultsSolver(system.getSolver()->solve(h)); break; // Simple simulation with numerical derivative
        case 2: functions.solve_varying_alpha(alpha_min, alpha_max, N_alpha, alpha_to_file); break;
        case 3: functions.solve_varying_dt(dt_min, dt_max, N_dt, dt_to_file); break;
        case 4: functions.solve_varying_N(Ns, N_to_file); break;
        case 5: best_alpha = functions.gradientDescent(initial_alpha, gamma, tolerance, Nmax_gradient, Nsteps_gradient); 
                cout << scientific << setprecision(5) << "best alpha= " << best_alpha << endl; break;
        case 6: system.getSolver()->thermalize();
                functions.printResultsSolver(system.getSolver()->solve(r_max, Nbins)); break;

        case 7: functions.runParallel(); break;
    } 
    
    auto stop = chrono::steady_clock::now(); // Store starting time to measure run time
    auto diff = stop - start; // Time difference
    cout << endl << fixed << "Simulation termined. Total running time: " << chrono::duration <double, milli> (diff).count()/1000 << " s" << endl; // Print run time*/
}