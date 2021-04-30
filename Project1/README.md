# Variational Monte Carlo for trapped bosons

To download the project go to the [main page of this repository](https://github.com/Matteo294/FYS4411) and click on the green button **Code** on the top right, then **Download ZIP**.  
Unzip the folder and open a terminal into folder **Project1**. Then cd into the folder **code**

```
cd code
```

Make the project by entering the command
```
make 

```
At this point the program can be run in multpliple ways, briefly listed here
0. **Single analytical simulation**: run a single simulation using analytical formulas for the local energy
1. **Single numerical simulation**: run a single simulation using the numerical approximation for the second derivative of the wavefunction when evaluating the local energy
2. **Varying dt**: run multiple simulations varying the time step length used to numerically integrate the Fokker-Planck equation
3. **Varying N**: run multiple simulation varying the number of particles
4. **Gradient descent**: run a gradient descent algorithm to find the best value for the variational parameter in the trial wavefunction. Then run a simple simulation using the obtained value for the variational parameter
5. **One body density**: simple simulation whose only goal is to store the data needed for the one-body density, skipping unrelevant calculations.

To run one the program in one particular mode one just needs to enter the corresponind id when launching the main. For example, in the case of the simple simulation (id 0) the command becomes
```
./main 0
```
All the relevant settings and parameters can be set in **main.cpp**.  
  
By default the flag **TO_FILE** is set to 0. This means that no data is saved to file and no further analysis can be carried on after the simulations is over (e.g. Blocking analysis). In order to store data to file we recommend to run first the bash script **setup.sh** which sets up all the necessary folders. In order to do this, run the commands
```
chmod +x setup.sh
sh setup.sh
```
At this point one can set the flag **TO_FILE** to 1 and proceed with the desired simulation. Some modalitites contain an extra flag for printing that can be set in the corresponding parameters section in **main.cpp**.  
  
