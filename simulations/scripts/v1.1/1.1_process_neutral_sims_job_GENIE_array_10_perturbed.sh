#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J dating_error_genie
#SBATCH -o logs/dating_error.%A_%a.out
#SBATCH -e logs/dating_error.%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH --mem=16G
#SBATCH -a 1-300

VERSION=1.1
DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}/replicates
GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

RMSE_VALUES=(0.27 0.30 0.33)
N_REPS=100

# task index -> (rmse, rep): 300 tasks = 3 rmse values x 100 reps
IDX=$((SLURM_ARRAY_TASK_ID - 1))
RMSE_IDX=$(( IDX / N_REPS ))
REP=$(( (IDX % N_REPS) + 1 ))
RMSE=${RMSE_VALUES[$RMSE_IDX]}

REP_DIR=${DATA}/rep${REP}
ANNOT=${REP_DIR}/dating_error/1.1_annotations_age_bins.rmse${RMSE}.txt
OUT=${REP_DIR}/dating_error/1.1_out_GENIE.rmse${RMSE}

echo "REP=${REP} RMSE=${RMSE}"

micromamba run -n GENIE ${GENIE} \
    --genotype ${REP_DIR}/1.1_neutral_out \
    --phenotype ${REP_DIR}/1.1_phenotypes.GENIE.txt \
    --annot ${ANNOT} \
    --output ${OUT} \
    --model G \
    --num-jack 20 \
    --verbose 1