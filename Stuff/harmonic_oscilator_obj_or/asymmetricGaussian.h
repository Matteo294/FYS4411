#include "wavefunction.h"
#include "math.h"

class AsymmetricGaussian: public Wavefunction{
    public:
        AsymmetricGaussian(class System* s, double alpha, double beta); // Constructor
        virtual double evaluateAll();
        virtual double evaluateSing(int part_idx);
        virtual double evaluateSecondDerivative();
        virtual double numericalSecondDerivative(int part_idx, int direction, double h);
    protected:
        double alpha; // Variational parameter
        double beta;
};