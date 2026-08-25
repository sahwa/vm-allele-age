#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.0_stabilising_selection
#SBATCH -o 2.0_stabilising_selection.%A_%a.out
#SBATCH -e 2.0_stabilising_selection.%A_%a.err
#SBATCH -p long
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -a 1-3

SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v2.0
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0

V_S=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${SCRIPTS}/2.0_stabilising_selection_VS_sweep.txt)

OUTFILE=${DATA}/2.0_stabilising_selection_VS_${V_S}

slim \
    -d MU=1.44e-8 \
    -d PI_TARGET=0.01 \
    -d V_S=${V_S} \
    -d END_TICK=50000 \
    -d "outfile='${OUTFILE}'" \
    ${SCRIPTS}/2.0_stabilising_selection.slim