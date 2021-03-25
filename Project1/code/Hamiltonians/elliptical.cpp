#include "elliptical.h"

Elliptical::~Elliptical(){};

Elliptical::Elliptical(System* system, double omegaXY, double omegaZ) : Hamiltonian(system, (int) 2){
    this->setParameter(0, omegaXY);
    this->setParameter(1, omegaZ);
}


double Elliptical::LocalEnergyAnalytic(){
    double alpha = this->system->getWavefunction()->getParameter(0);
    double beta = this->system->getWavefunction()->getParameter(1);
    double a = this->system->getWavefunction()->getParameter(2);
    double mass = this->system->getParticles()[0]->getMass();
    int Nparticles = this->system->getNParticles();
    double check=0.0;

    int i=0, j=0, m=0;
    double tmp=0.0, res_j = 0.0, res=0.0;
    vector<double> pos_i(this->system->getDimension(), 0.0);
    vector<double> sum_m(this->system->getDimension(), 0.0);
    vector<double> pos_im(this->system->getDimension(), 0.0);
    

    for(i=0; i<Nparticles; i++){
        pos_i = this->system->getParticles()[i]->getPosition();
        pos_i.back() *= beta; // I modify pos_i, the only place where I use pos_i is inside the j-cycle

        res_j = 0.0;
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
                
                tmp = this->system->cdot(sum_m, this->system->relative_position[i][j]);
                tmp += -4 * alpha * this->system->cdot(pos_i, this->system->relative_position[i][j]) + 2 + (a - 2 * this->system->relative_distance[i][j]) / (this->system->relative_distance[i][j] - a);
                tmp *= a / pow(this->system->relative_distance[i][j], 2) / (this->system->relative_distance[i][j] - a);
                res_j += tmp;

            }

        }

        res += res_j;

    }
    res *= - 0.5 / mass;
    res += -2 * pow(alpha, 2) * this->system->r2(pow(beta, 2)) / mass + 0.5 * mass * pow(this->params[0], 2) * this->system->r2(pow( this->params[1]/this->params[0], 2));
    res += alpha * (2 + beta) * Nparticles / mass;

    return res;
}



double Elliptical::LocalEnergyNumeric(double h){
    double alpha = this->system->getWavefunction()->getParameter(0);
    double mass = this->system->getParticles()[0]->getMass(); 
    int i=0, j=0;
    double res = 0.0;
    
    for(i=0; i<this->system->getNParticles(); i++){
        for(j=0; j<this->system->getDimension(); j++){
            res += this->system->getWavefunction()->numericalSecondDerivative(i, j, h);
        }
    }

    res *= -0.5 / mass / this->system->getWavefunction()->evaluateAll();
    res += 0.5 * mass * pow(this->params[0], 2) * this->system->r2(pow(this->params[1]/this->params[0], 2));
    return res;
    
}