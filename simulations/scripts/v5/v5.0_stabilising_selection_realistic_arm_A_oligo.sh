#!/bin/bash
# ARM A — oligogenic, realistic selection strength.
#   PI_TARGET=0.002 -> V_M/V_P=1e-3, V_S/V_P=19.5, theta=1,280/rep, ~1.4h/job
#   Run this arm first: short jobs, no resume, no checkpoint pileup.
#
#   sbatch --array=1-600 run_armA_oligo.sh
#SBATCH -A visscher-wray.prj
#SBATCH --job-name=stabA
#SBATCH --mem=48G
#SBATCH --partition=short
#SBATCH --output=armA_%A_%a.out
#SBATCH --error=armA_%A_%a.err
#SBATCH --array=1-600

BASE=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations
CONFIG=${BASE}/scripts/v5/armA_jobs.tsv
SLIM_SCRIPT=${BASE}/scripts/v5/v5.1_stabilising_selection_realistic.slim
OUT_DIR=${BASE}/data/v5.1/armA_oligo
mkdir -p ${OUT_DIR}

# fixed parameters -- INTEGER literals (25e3 / 5e8 / 75e3 are floats in Eidos)
MU=1.28e-8
NE=25000
PI_TARGET=0.002
L=500000000
END_TICK=75000
SD_EFFECT=0.1

# pick this task's line (skipping # comments)
LINE=$(grep -v '^#' "${CONFIG}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n "${LINE}" ]] || { echo "No line ${SLURM_ARRAY_TASK_ID} in ${CONFIG}" >&2; exit 1; }
read -r TAG V_S MEAN_EFFECT NEUTRAL <<< "${LINE}"

OUT="${OUT_DIR}/${TAG}"
echo "=== armA ${SLURM_ARRAY_TASK_ID}: ${TAG}  V_S=${V_S} mean=${MEAN_EFFECT} neutral=${NEUTRAL}"

if [[ -s "${OUT}.summary.tsv" && -s "${OUT}.trees" ]]; then
  echo "    already complete, skipping"; exit 0
fi

slim -s $((20261117 + SLURM_ARRAY_TASK_ID)) \
    -d "MU=${MU}" -d "NE=${NE}" -d "PI_TARGET=${PI_TARGET}" -d "L=${L}" \
    -d "V_S=${V_S}" -d "NEUTRAL=${NEUTRAL}" \
    -d "MEAN_EFFECT=${MEAN_EFFECT}" -d "SD_EFFECT=${SD_EFFECT}" \
    -d "END_TICK=${END_TICK}" -d "CKPT_EVERY=${END_TICK}" -d "LOG_EVERY=500" \
    -d "outfile='${OUT}'" \
    "${SLIM_SCRIPT}"

echo "=== done: ${TAG} ==="
