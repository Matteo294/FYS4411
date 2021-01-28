#pragma once
#include "system.h"
#include <vector>

class Particle{
    public:
        class System* system;
        std::vector<double> position;
        void setMass(double m);
        void setPosition(std::vector<double> new_pos);
        void move(std::vector<double> delta_pos);
        std::vector<double> getPosition();
        Particle(class System* system);
        double mass;
};