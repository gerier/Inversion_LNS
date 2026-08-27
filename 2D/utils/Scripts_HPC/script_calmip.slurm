#!/bin/bash

#SBATCH -J crosscorr_kernel
#SBATCH -N 10
#SBATCH -n 360
#SBATCH --ntasks-per-node=36
#SBATCH --ntasks-per-core=1
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
##SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR

export OMP_NUM_THREADS=1
dirname=${SLURM_JOBID}

SIMUFILE=$1

if [ -z "$SIMUFILE" ]; then
    echo "Erreur : aucun fichier de simulation fourni"
    exit 1
fi

mpirun -np 360  xinfrasound $SIMUFILE > output_${SLURM_JOBID}.log

~      
