#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.1_stabilising_selection
#SBATCH -o 2.1_stabilising_selection.%A_%a.out
#SBATCH -e 2.1_stabilising_selection.%A_%a.err
#SBATCH -p long
#SBATCH -c 8
#SBATCH --mem 128G
#SBATCH -a 1-4

VERSION="2.1"
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v${VERSION}
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}

V_S="$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${SCRIPTS}/${VERSION}_stabilising_selection_VS_sweep.txt)"
NE=25000
END_TICK=$(( NE * 15 ))

if [[ -z "${V_S}" ]]; then
    echo "ERROR: V_S empty for task ${SLURM_ARRAY_TASK_ID}"; exit 1
fi

echo "Running SLiM: NE=${NE}, V_S=${V_S}, L=5e8, END_TICK=${END_TICK}"
OUTFILE=${DATA}/${VERSION}_stabilising_selection_VS_${V_S}_NE_${NE}

stdbuf -oL slim \
    -d MU=1.44e-8 \
    -d PI_TARGET=0.01 \
    -d V_S=${V_S} \
    -d NE=${NE} \
    -d L=5e8 \
    -d END_TICK=${END_TICK} \
    -d "outfile='${OUTFILE}'" \
   s

