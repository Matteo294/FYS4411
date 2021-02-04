main : main.o system.o particle.o wavefunction.o gaussian.o asymmetricGaussian.o hamiltonian.o spherical.o solver.o metropolis.o importanceSampling.o random_generator.o
	g++ -Wall -o main main.o system.o particle.o wavefunction.o gaussian.o asymmetricGaussian.o hamiltonian.o spherical.o solver.o metropolis.o importanceSampling.o random_generator.o
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
asymmetricGaussian.o : asymmetricGaussian.h asymmetricGaussian.cpp
	g++ -c asymmetricGaussian.cpp
hamiltonian.o : hamiltonian.h hamiltonian.cpp
	g++ -c hamiltonian.cpp
spherical.o : spherical.h spherical.cpp
	g++ -c spherical.cpp
solver.o : solver.h solver.cpp
	g++ -c solver.cpp
metropolis.o : metropolis.h metropolis.cpp
	g++ -c metropolis.cpp
importanceSampling.o : importanceSampling.h importanceSampling.cpp
	g++ -c importanceSampling.cpp
random_generator.o : random_generator.h random_generator.cpp
	g++ -c random_generator.cpp
