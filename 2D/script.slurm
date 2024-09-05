#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --begin=now
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
#SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR
##  #SBATCH --partition=medium

module load openmpi/4.1.6-gcc12.2.0

mpirun -np 8 xinfrasound

