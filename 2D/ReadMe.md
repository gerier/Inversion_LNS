# Objective

This software was created for:
- infrasound wave propagation modelling
- sensitivity kernel of infrasound to background atmospheric model perturbations
- inversion of the background atmospheric model

The linearised Euler equations are at the heart of the three problems. To model them, finite differences method have been implemented. 

The file  `parameters.f90` contains a option to choose which application has to be executed.
In addition, this file contains all the parameters to design a simulation : time step, spatial step, domain specification, atmospheric model... 


# Execution

## General command
First, edit the file `parameters.f90`.

To compile the program:
```
make clean; make
```

To execute the program:
```
mpirun -np X inversion
```
with `X` the number of processus you may want to use. 

## Application on Pando (ISAE supercomputer)
First, edit the file `parameters.f90`.

Load specific packages:
```
module load openmpi/4.0.0-gcc8.2
```

To compile the program:
```
make clean; make
```

To execute the program:
Edit the file scrit.slurm with the number of processus, email, time of execute, etc.
```
sbatch script.slurm
```

## Application on Calmip supercomputer
First, edit the file `parameters.f90`.

Load specific packages:
```
module load gcc/7.3.0 openmpi/gnu/2.0.2.10
module unload gcc
```

To compile the program:
```
make clean; make
```

To execute the program:
Edit the file scrit.slurm with the number of processus, email, time of execute, etc.
```
sbatch script.slurm
```

# Organisation

## Files organisation

In the current directory, you may file to compile and execute the program: 
1. `parameters.f90`: file defining the different parameters [to be modified by the user]
2. `Makefile`: to compile easily the software
3. `xinfrasound`: the file to execute to run the software
3 directories contains the src files.

3 directories are present and contains the source files of the software [to not be modified by users]
1. `src`: files necessary to execute model the wave propagation
2. `src_gradient`: files necessary to execute the sensitivity kernel application (in addition to those in the folder `src`)
3. `src_inversion`: files necessary to execute the inversion application (in addition to those in the folder `src` and `src_gradient`)

In the directory `DATA`, you may find input files such as atmospheric models. 

3 directories contains examples of each application. 
TODO


## How the code works

![Texte alternatif](./doc/image_code.svg "Scheme of the software execution")







