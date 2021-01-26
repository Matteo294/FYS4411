main : main.o system.o particle.o wavefunction.o gaussian.o
	g++ -Wall -o main main.o system.o particle.o wavefunction.o gaussian.o
main.o : main.cpp
	g++ -c main.cpp
system.o : system.h system.cpp 
	g++ -c system.cpp
particle.o : particle.h particle.cpp
	g++ -c particle.cpp
wavefunction.o : wavefunction.h wavefunction.cpp
	g++ -c wavefunction.cpp
gaussian.o : gaussian.h gaussian.cpp
	g++ -c gaussian.cpp
