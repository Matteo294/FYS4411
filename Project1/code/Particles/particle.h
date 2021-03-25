#pragma once
#include "../System/system.h"
#include <vector>
using namespace std;

class Particle{
    public:
        Particle(class System* system, double mass, vector<double> pos);
        ~Particle();

        class System* system;

        /// Sets the mass of the particle
        void setMass(double m);
        /**Set the position of a particle**/
        /// Sets the position of the particle to the selected vector
        void setPosition(vector<double> new_pos);

        /// Returns a vector describing the position of the particle
        vector<double> getPosition();
        /// Returns the mass of the particle
        double getMass() {return this->mass;}
               
        /// This function varies the position of the particle of the vector delta_pos.
        void move(vector<double> delta_pos);

    private:
        double mass;
        vector<double> position;
};