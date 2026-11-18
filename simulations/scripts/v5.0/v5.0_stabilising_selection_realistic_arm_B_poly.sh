#!/bin/bash
# ARM B — polygenic, unrealistically strong selection. The SENSITIVITY CHECK.
#   PI_TARGET=0.10 -> V_M/V_P=1e-3, V_S/V_P=0.39, theta=64,000/rep, ~71h/job
#   lambda is identical to arm A, so the bias should be too.
#   Run only after arm A has produced a number.
#
#   sbatch --array=1-14 run_armB_poly.sh
#   requeue a timed-out task: sbatch --array=<id> run_armB_poly.sh
#SBATCH --job-name=stabB
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=long           # CHECK your cluster's partition name
#SBATCH --output=armB_%A_%a.out
#SBATCH --error=armB_%A_%a.err

set -euo pipefail

BASE=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations
CONFIG=${BASE}/scripts/v5.0/armA_jobs.tsv
SLIM_SCRIPT=${BASE}/scripts/v5.0/v5.0_stabilising_selection_realistic.slim
OUT_DIR=${BASE}/data/v5.0/armA_oligo


MU=1.28e-8
NE=25000
PI_TARGET=0.10
L=500000000
END_TICK=75000
SD_EFFECT=0.1
CKPT_EVERY=2500                    # ~2.4h of work per checkpoint

LINE=$(grep -v '^#' "${CONFIG}" | sed -n "${SLURM_ARRAY_TASK_ID}p")
[[ -n "${LINE}" ]] || { echo "No line ${SLURM_ARRAY_TASK_ID} in ${CONFIG}" >&2; exit 1; }
read -r TAG V_S MEAN_EFFECT NEUTRAL <<< "${LINE}"

OUT="${OUT_DIR}/${TAG}"
echo "=== armB ${SLURM_ARRAY_TASK_ID}: ${TAG}  V_S=${V_S} mean=${MEAN_EFFECT} neutral=${NEUTRAL}"

if [[ -s "${OUT}.summary.tsv" && -s "${OUT}.trees" ]]; then
  echo "    already complete, skipping"; exit 0
fi

# resume from the newest checkpoint, if one exists
RESUME_ARG=()
LATEST=$(ls -1 "${OUT}".tick*.trees 2>/dev/null | sort -V | tail -n1 || true)
if [[ -n "${LATEST}" ]]; then
  echo "    resuming from ${LATEST}"
  RESUME_ARG=(-d "RESUME='${LATEST}'")
fi

slim -s $((20261117 + 1000 + SLURM_ARRAY_TASK_ID)) \
    -d "MU=${MU}" -d "NE=${NE}" -d "PI_TARGET=${PI_TARGET}" -d "L=${L}" \
    -d "V_S=${V_S}" -d "NEUTRAL=${NEUTRAL}" \
    -d "MEAN_EFFECT=${MEAN_EFFECT}" -d "SD_EFFECT=${SD_EFFECT}" \
    -d "END_TICK=${END_TICK}" -d "CKPT_EVERY=${CKPT_EVERY}" -d "LOG_EVERY=1000" \
    "${RESUME_ARG[@]}" \
    -d "outfile='${OUT}'" \
    "${SLIM_SCRIPT}"

# keep only the newest checkpoint (15 x 590MB x 14 jobs would be ~120GB)
if [[ -s "${OUT}.trees" ]]; then
  rm -f "${OUT}".tick*.trees
else
  KEEP=$(ls -1 "${OUT}".tick*.trees 2>/dev/null | sort -V | tail -n1 || true)
  [[ -n "${KEEP}" ]] && ls -1 "${OUT}".tick*.trees | grep -v "^${KEEP}$" | xargs -r rm -f
fi

echo "=== done: ${TAG} ==="
