# Variational Monte Carlo for trapped bosons

To download the project go to the [main page of this repository](https://github.com/Matteo294/FYS4411) and click on the green button **Code** on the top right, then **Download ZIP**.  
Unzip the folder and open a terminal into folder **Project1**. Then cd into the folder **code**

```
cd code
```

**Make** the project by entering the command
```
make 
```
At this point the program can be run in multpliple ways, briefly listed here  
1. **Single analytical simulation** (ID 0): run a single simulation using analytical formulas for the local energy
2. **Single numerical simulation** (ID 1): run a single simulation using the numerical approximation for the second derivative of the wavefunction when evaluating the local energy
3. **Varying alpha** (ID 2): run multiple simulations varying the variational parameter in the wavefunction
4. **Varying dt** (ID 3): run multiple simulations varying the time step length used to numerically integrate the Fokker-Planck equation
5. **Varying N** (ID 4): run multiple simulations varying the number of particles
6. **Gradient descent** (ID 5): run a gradient descent algorithm to find the best value for the variational parameter in the trial wavefunction. Then run a simple simulation using the obtained value for the variational parameter
7. **One body density** (ID 6): simple simulation whose only purpose is to store the data needed for the one-body density, skipping unrelevant calculations.

To run the program in one particular mode one just needs to enter the corresponing ID when launching the main. For example, in the case of the simple simulation (ID 0) the command becomes
```
./main 0
```
All the relevant settings and parameters can be set in **main.cpp**. Read the documentation file **documentation.pdf** for detailed information.

By default the flag **RUN_PARALLEL** is set to 1. This means that the program is run in the parallelized version over the set number of threads (4 by defaul). The default number of threads can be modified by setting the parameter **NTHREADS**: we recommend to use a power of 2 as number of threads. Instead, by setting **RUN_PARALLEL** to 0, the program is run in the serial version.
  
By default the flag **TO_FILE** is set to 0. This means that no data is saved to file and no further analysis can be carried on after the simulations is over (e.g. Blocking analysis). In order to store data to file we recommend to run first the bash script **setup.sh** which sets up all the necessary folders. In order to do this, run the commands
```
chmod +x setup.sh
sh setup.sh
```
At this point one can set the flag **TO_FILE** to 1 and proceed with the desired simulation. Some modalitites contain an extra flag for printing to file that can be set in the corresponding parameters section in **main.cpp**.  
  
The Blocking variance analysis is carried by the script **blocking.py** in the **Analysis** folder. Once moved into the folder with the command
```
cd Analysis
```
it is possible to get a guide on how to run the script by entering
```
python blocking.py -h
```
or
```
python3 blocking.py -h
```
depending on your system
