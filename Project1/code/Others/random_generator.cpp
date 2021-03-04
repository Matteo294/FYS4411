#include "random_generator.h"

RandomGenerator::RandomGenerator(){

    uniform_real_distribution<double> uniform(0.0, 1.0);
    normal_distribution<double> normal(0.0, 1.0);
    this->uniform = uniform;
    this->normal = normal;
}