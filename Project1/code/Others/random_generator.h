#pragma once
#include <random>
using namespace std;

class RandomGenerator{
    
    public:
        RandomGenerator();
        uniform_real_distribution<double> uniform;
        normal_distribution<double> normal;
};
