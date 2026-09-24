#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.1_stabilising_selection
#SBATCH -o 2.1_stabilising_selection.%A_%a.out
#SBATCH -e 2.1_stabilising_selection.%A_%a.err
#SBATCH -p short
#SBATCH -c 4
#SBATCH -a 1-4

set -euo pipefail

VERSION="2.1"
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v${VERSION}
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SLIM_SCRIPT=${SCRIPTS}/${VERSION}_stabilising_selection.slim   # adjust to your filename

mkdir -p "${DATA}"

V_S="$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR == i' \
       "${SCRIPTS}/${VERSION}_stabilising_selection_VS_sweep.txt")"
if [[ -z "${V_S}" ]]; then
    echo "ERROR: V_S empty for task ${SLURM_ARRAY_TASK_ID}"; exit 1
fi

NE=25000
# END_TICK=$(( NE * 15 ))
# OUTFILE=${DATA}/${VERSION}_stabilising_selection_VS_${V_S}_NE_${NE}

# resume from the most recent checkpoint if one exists
# LATEST=$(ls -1v "${OUTFILE}".tick*.trees 2>/dev/null | tail -n 1) || true

# echo "NE=${NE}  V_S=${V_S}  L=5e8  END_TICK=${END_TICK}"
# echo "Resume from: ${LATEST:-<fresh start>}"

# stdbuf -oL slim \
#     -d MU=1.44e-8 \
#     -d PI_TARGET=0.01 \
#     -d V_S=${V_S} \
#     -d NE=${NE} \
#     -d L=5e8 \
#     -d END_TICK=${END_TICK} \
#     -d CKPT_EVERY=10000 \
#     -d "outfile='${OUTFILE}'" \
#     -d "RESUME='${LATEST:-}'" \
#     "${SLIM_SCRIPT}"

# Then we would like to recapitate and strip biallelics

# find the treee with the highest tick count

treefile=$(ls ${DATA}/2.1_stabilising_selection_VS_${VS}_NE_${NE}.tick*.trees | sort -V | tail -n1)
outfile=$(echo $treefile | sed 's/.trees/.biallelic.recapitated.trees/g')

micromamba run -n slim python 2.1_stabilising_selection_process_msprime_recapitate.py \
--tree ${treefile} \
--ancestral-ne ${NE} \
--rec-rate 1e-8 \
--seed ${SLURM_ARRAY_TASK_ID} \
--out ${outfile}
