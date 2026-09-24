#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.1_stabilising_selection
#SBATCH -o 2.1_stabilising_selection.%A_%a.out
#SBATCH -e 2.1_stabilising_selection.%A_%a.err
#SBATCH -p long
#SBATCH -c 1
#SBATCH --mem 128G
#SBATCH -a 1-4

set -euo pipefail

VERSION="2.0"
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v${VERSION}
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SLIM_SCRIPT=${SCRIPTS}/${VERSION}_stabilising_selection.slim   # adjust to your filename
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v${VERSION}

GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

REP_BASE=${DATA}/replicates/VS_${V_S}_NE_${N_E}/${SLURM_ARRAY_TASK_ID}

micromamba run -n slim python blah.py blah.py

micromamba run -n GENIE  ${GENIE} \
 --genotype ${REP_BASE}/${VERSION}_stabilising_selection_VS_${V_S}_NE_${N_E}.out \
 --phenotype ${REP_BASE}/${VERSION}_stabilising_selection_VS_${V_S}_NE_${N_E}_phenotypes.GENIE.txt \
 --annot ${REP_BASE}/${VERSION}_stabilising_selection_VS_${V_S}_NE_${N_E}_annotations_age_bins.txt \
 --output ${REP_BASE}/${VERSION}_stabilising_selection_VS_${V_S}_NE_${N_E}_neutral_out_GENIE \
 --model G \
 --verbose 1 \
 --nthreads 2


