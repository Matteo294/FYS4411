#ifndef PARTICLE_H
#define PARTICLE_H
/* Here is the definition and implementation of the class Particle. 
As the name suggests, the definition part simply consists in defining things without specifying what they make */
class Particle{
    private:
        double x; // Position
        double m; // Mass
        double v; // Speed
    public:
        Particle(double x, double m, double v); // Constructor: this function is run as soon as an object of the class is created
        void getInfo(); // This function will print information about the particle
        void accelerate();
};
#endif