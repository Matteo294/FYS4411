#include "elliptical.h"

Elliptical::~Elliptical(){};

Elliptical::Elliptical(System* system, double omegaXY, double omegaZ) : Hamiltonian(system){
    this->omegaXY = omegaXY;
    this->omegaZ = omegaZ;
}



double Elliptical::LocalEnergyAnalytic(){
    // things I need multiple times in the calculations or long expressions: better calculate them once for all
    // Hardcoded mass=1
    double alpha = this->system->getWavefunction()->getParameter(0);
    double beta = this->system->getWavefunction()->getParameter(1);
    double a = this->system->getWavefunction()->getParameter(2);
    int Nparticles = this->system->getNParticles();

    int i=0, j=0, m=0;
    double res = 0.0;
    vector<double> pos_i(this->system->getDimension(), 0.0);
    vector<double> sum_m(this->system->getDimension(), 0.0);
    vector<double> pos_im(this->system->getDimension(), 0.0);
    //vector<double> pos_ij(this->system->getDimension(), 0.0);
    

    for(i=0; i<Nparticles; i++){
        pos_i = this->system->getParticles()[i]->getPosition();
        pos_i.back() *= beta; // I modify pos_i, the only place where I use pos_i is inside the j-cycle

        for(j=0; j<Nparticles; j++){
            if(j!=i){
                               
                fill(sum_m.begin(), sum_m.end(), (double) 0.0); // reset sum_m to {0,0,0}

                for(m=0; m<Nparticles; m++){
                    if(m!=i){
                        pos_im = this->system->relative_position[i][m];
                        transform(pos_im.begin(), pos_im.end(), pos_im.begin(), bind1st(multiplies<double>(), a / pow(this->system->relative_distance[i][m],2) / (this->system->relative_distance[i][m] - a)));
                        transform(sum_m.begin(), sum_m.end(), pos_im.begin(), sum_m.begin(), plus<double>());
                    }
                }
                
                res += this->system->cdot(sum_m, this->system->relative_position[i][j]);
                res += -4 * alpha * this->system->cdot(pos_i, this->system->relative_position[i][j]) + 2 + (a - 2 * this->system->relative_distance[i][j]) / (this->system->relative_distance[i][j] - a);
                res *= a / pow(this->system->relative_distance[i][j], 2) / (this->system->relative_distance[i][j] - a);

            }

        }

        res *= - 0.5;

    }

    res += -2 * pow(alpha, 2) * this->system->r2(pow(beta, 2)) + 0.5 * this->system->r2(pow(this->getOmegaZ(), 2));
    res += alpha * (2 + beta) * Nparticles;

    return res;
}

double Elliptical::LocalEnergyNumeric(double h){return 0.0;}