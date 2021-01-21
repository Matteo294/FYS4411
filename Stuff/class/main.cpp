#include "particle.h"
#include "system.h"
#include <iostream>
#include <vector>
#include <time.h>

using namespace std;


double random(int min, int max){
    return ((double) rand()/RAND_MAX) * (max - min) + min; // random number between min and max
}

int main(){

    srand(time(NULL));

    /* ------------------ Program 1 ---------------------- */
    cout << "Part 1" << endl;

    /* We create an object of the class particle, assign its mass, position and velocity.
    Then we print those values */

    Particle p(0.0, 1.0, 1.5); // position 0, mass 1, velocity 1.5. As soon as we create the variable the constructor is called. p is one object of the class Particle.
    p.getInfo(); // call the function implemented in the Partice class that prints the attributes of the object (p)
    cout << endl << endl;
    /* --------------------------------------------------- */







    /* ------------------ Program 2 ---------------------- */
    cout << "Part 2" << endl;
    /* Here we make an array of class objects, assign parameters randomly and print them. 
    Then we accelerate and print the results */

    Particle p1(random(-5, 5), random(0, 10), random(-5, 5));
    Particle p2(random(-5, 5), random(0, 10), random(-5, 5));
    Particle p3(random(-5, 5), random(0, 10), random(-5, 5));
    p1.getInfo();
    p2.getInfo();
    p3.getInfo();
    cout << endl;

    // Vectors -> faster way to acces data
    vector<Particle> pvec{p1, p2, p3};
    // Printing becomes like this
    for(int i=0; i<3; i++){
        pvec[i].getInfo();
    }
    cout << endl;

    // Accelerate
    for(int i=0; i<3; i++){
        pvec[i].accelerate();
    }

    // Print results
    for(int i=0; i<3; i++){
        pvec[i].getInfo();
    }
    cout << endl << endl;  
    /* --------------------------------------------------- */









    /* ------------------ Program 3 ---------------------- */
    cout << "Part 3" << endl;
    /* Same as before but with a system class */
    Particle q1(random(-5, 5), random(0, 10), random(-5, 5));
    Particle q2(random(-5, 5), random(0, 10), random(-5, 5));
    Particle q3(random(-5, 5), random(0, 10), random(-5, 5));

    System s;

    s.addParticle(q1);
    s.addParticle(q2);
    s.addParticle(q3);

    s.getInfo();
    cout << endl;
    
    s.accelerateAll();
    
    s.getInfo();
    /* --------------------------------------------------- */



    return 0;
}