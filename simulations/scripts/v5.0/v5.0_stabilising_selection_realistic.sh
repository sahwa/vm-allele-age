#!/bin/bash
#SBATCH --job-name=stabsel
#SBATCH --array=1-30
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G                 # CHECK: scales with NE x segregating sites
#SBATCH --time=48:00:00           # CHECK against your partition limits
#SBATCH --partition=long          # CHECK: your cluster's partition name
#SBATCH --output=logs/stabsel_%A_%a.out
#SBATCH --error=logs/stabsel_%A_%a.err

set -euo pipefail

OUT_DIR=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v5.0
SLIM_SCRIPT=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/stabsel_v3.slim
mkdir -p "${OUT_DIR}" logs

# ---------------------------------------------------------------------------
# Fixed parameters. INTEGER LITERALS, NOT 25e3 / 5e8 / 75e3 -- those are floats
# in Eidos and addSubpop(), initializeGenomicElement() and the END_TICK tick
# specifier all want integers. This is the most likely way the original call
# fails, and it fails late.
# ---------------------------------------------------------------------------
MU=1.28e-8                        # per base per generation (float is fine here)
NE=25000
PI_TARGET=0.10
L=500000000
END_TICK=75000
SD_EFFECT=0.1

# ---------------------------------------------------------------------------
# Job table:  name  V_S  MEAN_EFFECT  NEUTRAL  REPS
#
# V_S is passed even for the neutral runs so that every constant the SLiM
# script references is defined (LAMBDA_PRED touches V_S).
#
# lambda_pred = sigma2_m / V_S, with sigma2_m = MEAN_EFFECT^2 + SD_EFFECT^2.
# Note the skew runs therefore have a LARGER sigma2_m, so their lambda is not
# the same as the symmetric run at the same V_S -- do not compare them as if
# lambda were held fixed.
# ---------------------------------------------------------------------------
CONFIGS=(
  "null         10    0.00   T   4"    # lambda = 0. MUST return bias = 1.000
  "realistic    10    0.00   F   8"    # lambda ~ 1e-3, the number you quote
  "loud          1    0.00   F   5"    # lambda ~ 1e-2, verify functional form
  "veryloud    0.1    0.00   F   5"    # lambda ~ 1e-1, find where theory breaks
  "skewmild     10    0.05   F   4"    # E[beta] != 0, ~20% directional
  "skewstrong   10    0.10   F   4"    # E[beta] != 0, ~50% directional
)

# expand configs x reps into a flat job list
JOBS=()
for cfg in "${CONFIGS[@]}"; do
  read -r name vs mean neut reps <<< "${cfg}"
  for ((r = 1; r <= reps; r++)); do
    JOBS+=("${name} ${vs} ${mean} ${neut} ${r}")
  done
done

N_JOBS=${#JOBS[@]}

# `sbatch --array` is set by hand above, so fail loudly if they disagree
if [[ "${SLURM_ARRAY_TASK_ID:-0}" -gt "${N_JOBS}" ]]; then
  echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} exceeds N_JOBS=${N_JOBS}" >&2
  exit 1
fi
echo "Job table holds ${N_JOBS} jobs; set --array=1-${N_JOBS}"

read -r NAME V_S MEAN_EFFECT NEUTRAL REP <<< "${JOBS[$((SLURM_ARRAY_TASK_ID - 1))]}"

TAG="${NAME}_rep${REP}"
OUT="${OUT_DIR}/${TAG}"

echo "=== task ${SLURM_ARRAY_TASK_ID}/${N_JOBS}: ${TAG} ==="
echo "    V_S=${V_S}  MEAN_EFFECT=${MEAN_EFFECT}  NEUTRAL=${NEUTRAL}"

# ---------------------------------------------------------------------------
# Resume from the newest checkpoint if one exists. END_TICK=75000 at L=5e8 may
# not finish inside a single wall-clock allocation; requeue the same array task
# and it picks up where it stopped.
# ---------------------------------------------------------------------------
RESUME_ARG=()
LATEST=$(ls -1 "${OUT}".tick*.trees 2>/dev/null | sort -t. -k2 -V | tail -n1 || true)
if [[ -n "${LATEST}" ]]; then
  echo "    resuming from ${LATEST}"
  RESUME_ARG=(-d "RESUME='${LATEST}'")
fi

# different seed per task, reproducible from the array id
SEED=$((20261117 + SLURM_ARRAY_TASK_ID))

slim \
    -s "${SEED}" \
    -d "MU=${MU}" \
    -d "NE=${NE}" \
    -d "PI_TARGET=${PI_TARGET}" \
    -d "L=${L}" \
    -d "V_S=${V_S}" \
    -d "NEUTRAL=${NEUTRAL}" \
    -d "MEAN_EFFECT=${MEAN_EFFECT}" \
    -d "SD_EFFECT=${SD_EFFECT}" \
    -d "END_TICK=${END_TICK}" \
    -d "CKPT_EVERY=5000" \
    -d "LOG_EVERY=1000" \
    "${RESUME_ARG[@]}" \
    -d "outfile='${OUT}'" \
    "${SLIM_SCRIPT}"

echo "=== done: ${TAG} ==="

# collect afterwards with:
#   head -n1 ${OUT_DIR}/null_rep1.summary.tsv > ${OUT_DIR}/all_summaries.tsv
#   tail -qn1 ${OUT_DIR}/*.summary.tsv       >> ${OUT_DIR}/all_summaries.tsv
