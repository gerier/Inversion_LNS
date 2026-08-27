#!/bin/bash

#SBATCH -J script_test_roland
#SBATCH --nodes=16
#SBATCH --ntasks=384
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=solene.gerier@isae-supaero.fr
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=test
##SBATCH -o oslurm.%j.out # STDOUT
#SBATCH -e eslurm.%j.err # STDERR
#SBATCH --partition=medium
dirname=${SLURM_JOBID}

SIMUFILE=$1

if [ -z "$SIMUFILE" ]; then
    echo "Erreur : aucun fichier de simulation fourni"
    exit 1
fi

mpirun -np 384 xinfrasound $SIMUFILE > output_${SLURM_JOBID}.log

~      
