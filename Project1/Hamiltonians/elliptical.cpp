#include "elliptical.h"

Elliptical::Elliptical(System* system, double omegaXY, double omegaZ) : Hamiltonian(system){
    this->omegaXY = omegaXY;
    this->omegaZ = omegaZ;
}

double Elliptical::potential(){
    // Hardcoded mass=1
    double V = 0.0;
    for(int i=0; i<this->system->getNParticles(); i++){
        for(int j=0; j<2; j++){
            V += pow(this->getOmegaXY()*this->system->getParticles()[i]->getPosition()[j], 2);
        }
        V += pow(this->getOmegaZ()*this->system->getParticles()[i]->getPosition()[2], 2);
    }
    return 0.5*V;
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
    vector<double> pos_im(this->system->getDimension(), 0.0);
    vector<double> sum_m(this->system->getDimension(), 0.0);
    vector<double> pos_ij(this->system->getDimension(), 0.0);
    double dist_im=0.0;
    double dist_ij=0.0;


    for(i=0; i<Nparticles; i++){
        pos_i = this->system->getParticles()[i]->getPosition();
        pos_i.back() *= beta; // I modify pos_i, the only place where I use pos_i is inside the j-cycle

        for(j=0; j<Nparticles; j++){
            if(j!=i){
                pos_ij = this->system->getParticles()[i]->getRelativePosition(j);
                dist_ij = this->system->getParticles()[i]->getRelativeDistance(j);

                for(m=0; m<Nparticles; m++){
                    fill(sum_m.begin(), sum_m.end(), (double) 0.0); // reset sum_m to {0,0,0}
                    if(m!=i){
                        pos_im = this->system->getParticles()[i]->getRelativePosition(m);
                        dist_im = this->system->getParticles()[i]->getRelativeDistance(m);
                        transform(pos_im.begin(), pos_im.end(), pos_im.begin(), bind1st(multiplies<double>(), a / pow(dist_im,2) / (dist_im - a)));
                        transform(sum_m.begin(), sum_m.end(), pos_im.begin(), sum_m.begin(), plus<double>());
                    }
                }
                
                res += this->system->cdot(sum_m, pos_ij);
                res += -4 * alpha * this->system->cdot(pos_i, pos_ij) + 2 + (a - 2 * dist_ij) / (dist_ij - a);
                res *= a / pow(dist_ij, 2) / (dist_ij - a);

            }

        }

    }

    res += 4 * pow(alpha, 2) * this->system->r2(pow(beta, 2)) + this->system->r2(pow(this->getOmegaZ(), 2));

    return res/2;
}

double Elliptical::LocalEnergyNumeric(double h){return 0.0;}
