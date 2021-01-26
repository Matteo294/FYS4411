#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta); // Constructor
        double evaluate();
        double evaluateSecondDerivative();
        double numericalSecondDerivative();
    private:
        double alpha; // Variational parameter
        double beta;
};