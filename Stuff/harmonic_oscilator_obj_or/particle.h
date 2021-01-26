#pragma once
#include "system.h"
#include <vector>

class Particle{
    public:
        class System* system;
        std::vector<double> position;
        double mass;

        void setMass(double m);
        void setPosition(std::vector<double> new_pos);
        void updPosition(std::vector<double> var_pos);
        std::vector<double> getPosition();
        Particle();
};