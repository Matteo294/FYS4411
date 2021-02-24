#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta); // Constructor
        double evaluateAll();
        double evaluateSing(int part_idx);
        double analyticalSecondDerivative();
        double analyticalAlphaDerivative();
        double numericalSecondDerivative(int part_idx, int direction, double h);
        vector<double> DriftForce(int part_idx);
    protected:
        double alpha; // Variational parameter
        double beta; // Asymmetric term
};