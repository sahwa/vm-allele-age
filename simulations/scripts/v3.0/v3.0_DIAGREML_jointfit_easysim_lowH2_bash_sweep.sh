#!/bin/bash
#SBATCH --job-name=vm_v5
#SBATCH --partition=short
#SBATCH --array=1-300
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=v5_%A_%a.out
#SBATCH --error=v5_%A_%a.err

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

BASE=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations
Rscript v3.0_DIAGREML_jointfit_easysim_lowH2_bash_sweep.R cell $(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/data/v5.1/cells.txt)

