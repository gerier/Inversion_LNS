#!/bin/bash

#SBATCH --nodes=15
#SBATCH --ntasks=360
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --begin=now
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
#SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR
#SBATCH --partition=medium

module load python/3.8
# module load mpi4py/3.0.1-openmpi4.0.0-gcc8.2-python3.7
source activate env_df1d

mpirun -np 360 python3 ./solve_inverseproblem_2sources.py
