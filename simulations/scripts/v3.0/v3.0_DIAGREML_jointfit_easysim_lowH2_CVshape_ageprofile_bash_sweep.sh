#!/bin/bash
#SBATCH --job-name=DIAGREML_jointfit_easysim_lowH2_CVshape_bash_sweep
#SBATCH --partition=short
#SBATCH --array=1-25
#SBATCH --cpus-per-task=2
#SBATCH --output=DIAGREML_jointfit_easysim_lowH2_CVshape_bash_sweep_%A_%a.out
#SBATCH --error=DIAGREML_jointfit_easysim_lowH2_CVshape_bash_sweep_%A_%a.err

set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

BASE=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations
Rscript v3.0_DIAGREML_jointfit_easysim_lowH2_CVshape_ageprofile_bash_sweep.R cell $(sed -n "${SLURM_ARRAY_TASK_ID}p" /well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v3.4/cells.txt)

