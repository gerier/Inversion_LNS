#!/bin/bash


#SBATCH -J script_test_roland
#SBATCH -N 10
#SBATCH -n 100
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
#SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR

export OMP_NUM_THREADS=1
dirname=${SLURM_JOBID}

mpirun -np 100 xinfrasound  > output_${SLURM_JOBID}.log

~      
