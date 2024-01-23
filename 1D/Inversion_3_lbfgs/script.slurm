#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --begin=now
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
#SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR

module load python/3.8
# module load mpi4py/3.0.1-openmpi4.0.0-gcc8.2-python3.7 
source activate env_df1d

python3 ./solve_inverseproblem.py

