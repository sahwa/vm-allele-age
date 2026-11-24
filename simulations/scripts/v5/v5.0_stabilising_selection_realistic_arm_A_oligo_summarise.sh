#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH --job-name=stabA_summary
#SBATCH --mem=48G
#SBATCH --partition=short
#SBATCH --output=armA_summary_%A_%a.out
#SBATCH --error=armA_summary_%A_%a.err

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v5.0/armA_oligo 
OUT=${DATA}/v5.0_stabilising_selection_realistic_arm_A_oligo_summary.tsv

python3 v5.0_stabilising_selection_realistic_arm_A_oligo.py ${DATA} --out ${OUT}
