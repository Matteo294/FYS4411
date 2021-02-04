#pragma once
#include "system.h"
#include <vector>
using namespace std;

class Particle{
    public:

        // Constructor
        Particle(class System* system, double mass, vector<double> pos);
        
        // Setters
        void setMass(double m);
        void setPosition(vector<double> new_pos);

        // Getters
        vector<double> getPosition();
        double getMass() {return this->mass;}

        // Other functions
        void move(vector<double> delta_pos);

        // Attributes
        class System* system;

    private:
        // Private attributes: access only via setters and getters
        double mass;
        vector<double> position;
};