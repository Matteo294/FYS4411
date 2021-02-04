#include<iostream>
#include "particle.h"

/* Definitions
class = set of rules, functions, variables that are shared within different objects. We want to define this functions and parameters once for all and then
        create as many objects as we need that will automatically inherit this properties.
object = a variable that contains all the things defined in a class */


/* Here is the constructor implementation: this function is run as soon as an object is created.
In this case we simply want to store the parameters x, m, v passed to the function. In fact, parameters passed to
the function, are volatile and get destroyed at the end of the function. 

(Not important) By writing "this->x" we mean that we want
to acces the variable x inside this class. The syntax is an abbreviation for (*this).x, where this is an automatically
generated pointer to the object that will be created from this class. */

Particle::Particle(double x, double m, double v){
    this->x = x;
    this->m = m;
    this->v = v;
}


/* Here is the implementation of the function getInfo */

void Particle::getInfo(){
    std::cout << "Position: " << this->x << "  Velocity: " << this->v << "  Mass: " << this->m << std::endl;
}

/* Here is the implementation of the accelerate function that doubles the velocity */
void Particle::accelerate(){
    this->v *= 2;
}
/* ----------------------------------------------------------- */