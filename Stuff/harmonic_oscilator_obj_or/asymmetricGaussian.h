#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta); // Constructor
        virtual double evaluate();
        virtual double evaluateSecondDerivative();
        virtual double numericalSecondDerivative();
    protected:
        double alpha; // Variational parameter
        double beta;
};