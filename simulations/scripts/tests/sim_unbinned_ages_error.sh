#!/bin/bash
#SBATCH --job-name=vm_unbinned
#SBATCH --array=1-500
#SBATCH --cpus-per-task=1
#SBATCH --output=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/sim_unbinned_ages_err/logs/%A_%a.out
#SBATCH --error=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/sim_unbinned_ages_err/logs/%A_%a.err

## Usage:
##   mkdir -p .../data/tests/sim_unbinned_ages_err/logs
##   Rscript sim_unbinned_ages_error.R grid     # prints the cell count
##   sbatch --array=1-<cell count> sim_unbinned_submit.sh    # override the default above
##   Rscript scripts/tests/sim_ubinned_ages_error.R collect
##
## Memory: C is N_inds x n_t doubles inside fit_sel_reml. With N_inds = 1e5 and
## ~60 age bins that is ~50 MB, so 4G is generous; raise it if N_AGE_BINS or
## N_inds go up. Time: the error cells are the slow ones, since the kernel
## support widens with age.

set -euo pipefail

SCRIPT=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/tests/sim_unbinned_ages_error.R

Rscript "$SCRIPT" cell "${SLURM_ARRAY_TASK_ID}"
