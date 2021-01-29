#pragma once
#include "system.h"
#include <vector>

class Particle{
    public:

        // Constructor
        Particle(class System* system, double mass, std::vector<double> pos);

        // Attributes
        class System* system;
        
        // Setters
        void setMass(double m);
        void setPosition(std::vector<double> new_pos);

        // Getters
        std::vector<double> getPosition();
        double getMass() {return this->mass;}

        // Other functions
        std::vector<double> position;
        void move(std::vector<double> delta_pos);

    private:
        // Private attributes: access only via setters and getters
        double mass;
};