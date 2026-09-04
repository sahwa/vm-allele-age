#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.0_stabilising_selection
#SBATCH -o 2.0_stabilising_selection.%A_%a.out
#SBATCH -e 2.0_stabilising_selection.%A_%a.err
#SBATCH -p short
#SBATCH -c 10
#SBATCH -a 1

SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v2.0
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v2.0

# V_S=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${SCRIPTS}/2.0_stabilising_selection_VS_sweep.txt)
NE=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${SCRIPTS}/2.0_stabilising_selection_NE_sweep.txt)
V_S=5

echo "Running Slim using ${NE} individuals and setting V_S as ${V_S}"

OUTFILE=${DATA}/2.0_stabilising_selection_VS_${V_S}_NE_${NE}

#stdbuf -oL slim \
#    -d MU=1.44e-8 \
#    -d PI_TARGET=0.01 \
#    -d V_S=${V_S} \
#    -d NE=${NE} \
#    -d END_TICK=150000 \
#    -d "outfile='${OUTFILE}'" \
#    ${SCRIPTS}/2.0_stabilising_selection.slim


stdbuf -oL slim \
	-d MU=1.44e-8 \
	-d PI_TARGET=0.01 \
	-d V_S=5 \
	-d NE=25000 \
  	-d L=5e8 \
	-d END_TICK=2000 \
	-d "outfile='/tmp/timing_test.trees'" \
  	2.0_stabilising_selection.slim
