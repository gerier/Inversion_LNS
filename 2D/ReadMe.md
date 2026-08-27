# Objective

This software was developed for:
- infrasound wave propagation modelling
- computation of infrasound sensitivity kernel with respect to background atmospheric parameters
- inversion of the background atmospheric model

The linearised Euler equations are at the heart of the infrasound modelling. These equations are solved using the finite-difference method.

The `.parfile` file contains an option to select which application should be executed.
In addition, it contains all the parameters required to define a simulation, such as time step, spatial step, domain specification, atmospheric model, etc.


# Execution

## General command

To compile the program:
```
make clean; make
```

To execute the program:
```
mpirun -np X xinfrasound
```
where `X` the number of processes to use.

## Application on Pando (ISAE supercomputer)

First, edit the `example.parfile` file.

Load the required packages:
```
module load openmpi/4.0.0-gcc8.2
```

To execute the program, edit the `scrit_pando.slurm` file and specify the required parameters, sucb as the number of processes, email address, maximum execution time, etc.
Then run:
```
sbatch script_pando.slurm example.parfile
```

## Application on CALMIP supercomputer
First, edit the `example.parfile` file.

Load the required packages:
```
module load gcc/7.3.0 openmpi/gnu/2.0.2.10
module unload gcc
```
For more information, see the CALMIP  documentation at https://www.calmip.univ-toulouse.fr/espace-utilisateurs/doc-technique-olympe/ameliorer-les-performances/compilateurs-et-wrappers-mpi

To execute the program, edit the `scrit_calmip.slurm` file and specify the required parameters, such as the number of processes, email address, the maximum execution time, etc.
Then run:
```
sbatch script_calmip.slurm example.parfile
```

# Organisation

## Files organisation

The current directory contains the files required to compile and execute the program:
1. `example.parfile`: file defining the user parameters
2. `Makefile`: file used to compile the software
3. `xinfrasound`: executable used to run the software



The source code is organised into four directories. These files should not be modified by users.
1. `src`: main source files, including MPI parallelisation, output writing, initialisation, etc.
2. `src_forward`: files required to model wave propagation
3. `src_gradient`: files required to compute the sensitivity kernels, in addition to the files in the `src` directory
4. `src_inversion`: files required to perform the inversion, in addition to the files in the `src` and `src_gradient` directories

The `DATA` directory contains input files, such as atmospheric models or infrasound observations (pressure time series).

The `OUTPUT`directory contains the output of the run.

The `Example` directory contains the application cases present in Gerier et al. Sensitivity kernels of infrasound arrival-time to atmospheric parameters from adjoint method (in prep.).

# Define a simulation

The `.parfile` defines the following required parameters:


| Parameter | English description | Conditions / notes |
|---|---|---|
| `method` | Selects the type of computation: 1 = forward problem, 2 = sensitivity-kernel computation, 3 = inversion. | |
| `validation` | Controls how the sensitivity kernel is computed when `method = 2`: `.true.` uses small perturbations; `.false.` uses the adjoint method. | Used only for kernel computation. |
| `which_kernel` | Selects the sensitivity kernel to be computed. |  |
| `NPROC_X` | Number of MPI processes along the horizontal (Ox) direction. |  |
| `NPROC_Y` | Number of MPI processes along the vertical (Oy) direction. |  |
| `NX` | Number of grid points along the horizontal (Ox) direction. | Must be a multiple of `NPROC_X`. |
| `NY` | Number of grid points along the vertical (Oy) direction. | Must be a multiple of `NPROC_Y`. |
| `DELTAX` | Spatial grid spacing along (Ox), used by the forward problem. | Must satisfy the spatial dispersion criterion. |
| `DELTAY` | Spatial grid spacing along (Oy), used by the forward problem. | Must satisfy the spatial dispersion criterion. |
| `NSTEP` |  Number of time iterations performed for the forward problem. |  |
| `DELTAT` | Time step used for the forward problem. | Must satisfy the CFL stability criterion. |
| `USE_PML_XMIN` | Enables the Perfectly Matched Layer (PML) on the left boundary. | If both X-side PML flags are false, periodic boundary conditions are used on the left and right. |
| `USE_PML_XMAX` | Enables the Perfectly Matched Layer (PML) on the right boundary. | If both X-side PML flags are false, periodic boundary conditions are used on the left and right. |
| `USE_PML_YMIN` | Enables the PML on the lower boundary. | If both Y-side PML flags are false, periodic boundary conditions are used on the lower and upper boundaries. The parfile comment states that a false YMIN flag makes the lower boundary reflective. |
| `USE_PML_YMAX` | Enables the PML on the upper boundary. | If both Y-side PML flags are false, periodic boundary conditions are used on the lower and upper boundaries. |
| `f0` | Fundamental frequency of the source. |  |
| `type_source` | Selects the spatial source type: 1 = plane wave, 2 = point source, 3 = point source with SSF spreading. | Only used for a plane wave. |
| `wavefront` | Selects the propagation direction of a plane-wave source: 1 = x direction, 2 = y direction. | Used only when `type_source = 1`. |
| `xsource` | Horizontal (Ox) coordinate of the source. |  |
| `ysource` | Vertical (Oy) coordinate of the source. |  |
| `SSF_Sigma` | Spatial spreading parameter of the source. | Used only when `type_source = 3`; the parfile indicates that it is typically of the order of one wavelength. |
| `observation` | Selects the observation type used for kernel computation: 0 = full waveform, 1 = arrival time. |  |
| `observation_from_file` | Selects how observations are obtained: 0 = simulated, 1 = read from file, 2 = travel time used as a time shift when `observation = 2`, 3 = arbitrary time shift equal to 1. | The meanings are taken directly from the parfile comments. |
| `path_obs_file` | Path/prefix of the observation file. | Used when `observation_from_file = 1`. |
| `window_waveform` | Selects a time window (1) or take all the waveform (0) to compute the adjoint source. | Used when `observation = 0`. |
| `NREC_SET` | Number of receiver sets to define. One receiver set corresponds to one receiver line. |  |
| `NREC_PER_SET` | Number of receivers in each receiver set. | Array with `NREC_SET` entries. |
| `REC_SET_INFO_1` | Defines the start and end coordinates of receiver set 1: x_start, y_start, x_end, y_end. | One line is required for each receiver set. |
| `REC_wr_1` | Defines the receiver cross-correlation window: start time and duration. | One entry is specified per receiver. |
| `atmospheric_model_file` | Selects whether the atmospheric model is read from a file (`.true.`) or defined by hard-coded physical values (`.false.`). |  |
| `atmospheric_file_name_true` | Name/path of the true atmospheric model file. | Used when `atmospheric_model_file = .true.`; for the forward problem, this is the model that is used. |
| `atmospheric_file_name_prior` | Name/path of the prior atmospheric model file. | Used when `atmospheric_model_file = .true.`; the original documentation states that it is not used for `method = 1`. |
| `gamma_chemestry_value` | Ratio of specific heats used for the true/prior atmospheric model when no model file is used. | Used when `atmospheric_model_file = .false.`. The spelling follows the parameter name in the parfile. |
| `cp_unrelaxed_true` | Sound propagation speed of the true model when no atmospheric model file is used. | Used when `atmospheric_model_file = .false.`. |
| `density_true` | Density of the true model when no atmospheric model file is used. | Used when `atmospheric_model_file = .false.`. |
| `windx_value_true` |Wind speed along the x direction in the true model when no atmospheric model file is used. | Used when `atmospheric_model_file = .false.`. |
| `cp_unrelaxed_prior` | Sound propagation speed of the prior model when no atmospheric model file is used. | Used when `atmospheric_model_file = .false.`. |
| `density_prior` | Density of the prior model when no atmospheric model file is used. | Used when `atmospheric_model_file = .false.`. |
| `windx_value_prior` | Wind speed along the x direction in the prior model when no atmospheric model file is used. | Used when `atmospheric_model_file = .false.`. |
| `NPERTURB_MODEL` | Number of perturbations added to the atmospheric model. | 0 means that no perturbation is added. |
| `add_windperturb_profile` | Controls whether a Gaussian perturbation is added to the wind profile. | `.true.` enables the perturbation. |
| `ymin_wind` | Minimum altitude of the added Gaussian wind perturbation. | Used when `add_windperturb_profile = .true.`. |
| `ymax_wind` | Maximum altitude of the added Gaussian wind perturbation. | Used when `add_windperturb_profile = .true.`. |
| `mean_gauss_wind` | Mean altitude/position parameter of the Gaussian wind perturbation. | Used when `add_windperturb_profile = .true.`; the parfile comment gives km as the unit. |
| `sigma2_gauss_wind` | Variance parameter of the Gaussian wind perturbation. | Used when `add_windperturb_profile = .true.`; the parfile comment gives km² as the unit. |
| `max_wind_factor` | Strength (maximum amplitude) of the added Gaussian wind perturbation. | Used when `add_windperturb_profile = .true.`; the parfile comment gives m/s as the unit. |
| `IT_DISPLAY` | Frequency at which forward-problem images are saved and information about the forward simulation is displayed. |  |
| `save_normimage_overtime` | Controls normalization of saved forward-problem images: 0 = normalize with respect to the current iteration, 1 = normalize with respect to all iterations. | Used for `method = 1`. |
| `save_adjoint_source` | Controls whether the adjoint source is saved. |  |
| `N_GLOB_FRAMES` | Number of global checkpoint frames. |  |
| `N_LOC_FRAMES` | Number of local checkpoint frames. |  |
| `parametrisation` | Selects the inversion parameterization. The documented options include 1: rho+c+windx; 2: rho+p+windx; 3: log(rho)+log(c)+windx; 35: log(rho)+c+windx; 36: log(rho/rho_prior)+c+windx; 4: log(rho)+log(c)+windx; 5: log(c)+log(p)+windx. | Used for inversion (`method = 3`). |
| `scale_model` | Scaling factors applied to the three inversion variables. | Three values, in the order defined by `parametrisation`. |
| `type_gradient` | Selects the optimization descent direction: 1 = Fletcher-Reeves, 2 = Polak-Ribiere, 3 = Perry-Shanno, 4 = Hager-Zhang, 5 = Dai-Kou, 6 = L-BFGS. |  |
| `mem_lbfgs` | Memory size of the L-BFGS method. | Used when `type_gradient = 6`. |
| `type_regul_term` | Selects whether a regularization term is added: 0 = none, 1 = norm(m - m_prior). |  |
| `regul_weight` | Weight applied to the regularization term. |  |
| `type_smoothing` | Selects the smoothing applied after a model update: 1 = mean, 2 = Gaussian, 3 = median. |  |
| `steepest_nbiter_default` | Number of iterations in the first part of the optimization algorithm, which uses backtracking combined with steepest descent. |  |
| `maxiter_innerloop` | Maximum number of iterations used to search for the descent step length. |  |
| `maxiter_outerloop` | Maximum number of iterations of the optimization algorithm, i.e. the maximum number of models tested. |  |
| `tol_x` | Tolerance criterion for stopping the optimization when the model is no longer improved. | The original documentation states that this parameter is currently not used. |
| `alpha_max` | Maximum step length that can be tested during the step-length search. |  |
| `rate` | Step reduction factor used by the backtracking algorithm. |  |
| `c1` |  Parameter of the strong Wolfe conditions controlling whether the decrease in the objective function is sufficient. |  |
| `c2` | Parameter of the strong Wolfe conditions controlling the curvature condition. |  |
