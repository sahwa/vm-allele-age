#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 4.0_stabilising_selection
#SBATCH -o 4.0_stabilising_selection.%A_%a.out
#SBATCH -e 4.0_stabilising_selection.%A_%a.err
#SBATCH -p short
#SBATCH -a 1-100
#SBATCH --mem 8G

# array indices 0-99
KS=(0 0.002 0.005 0.01)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
K=${KS[$((IDX % 4))]}
SEED=$((IDX / 4 + 1))

OUT=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0

slim -s $SEED -d N=50000 -d N_SAMPLE=2500 -d PI_TARGET=1.0 -d END_TICK=600 \
     -d K_SEL=$K -d "OUTFILE_TICKS='/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0/sing_K${K}_seed${SEED}.TICKS.tsv'" -d "OUTFILE_END='/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v4.0/sing_K${K}_seed${SEED}.END.tsv'" v4.0_selection.slim
