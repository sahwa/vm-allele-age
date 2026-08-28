#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.0_stabilising_selection
#SBATCH -o 2.0_stabilising_selection.%A_%a.out
#SBATCH -e 2.0_stabilising_selection.%A_%a.err
#SBATCH -p long
#SBATCH -c 1
#SBATCH -a 1-10

SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v2.0
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0

micromamba run -n slim python ${SCRIPTS}/2.0_stabilising_selection_process.py ${SLURM_ARRAY_TASK_ID}