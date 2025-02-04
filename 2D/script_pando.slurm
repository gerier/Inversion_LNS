#!/bin/bash


#SBATCH -J script_test_roland
#SBATCH --nodes=6
#SBATCH --ntasks=144
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
##SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR
#SBATCH --partition=medium
dirname=${SLURM_JOBID}

mpirun -np 144 xinfrasound  > output_${SLURM_JOBID}.log

~      
