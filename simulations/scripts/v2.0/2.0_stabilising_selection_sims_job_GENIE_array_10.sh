#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.0_process_array
#SBATCH -o logs/2.0_process.%A_%a.out
#SBATCH -e logs/2.0_process.%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH -a 1-10

VERSION=2.0
SELECTION_TYPE="stabilising_selection"
VS=5
NE=20000

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v${VERSION}

GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

micromamba run -n slim python ${SCRIPTS}/2.0_stabilising_selection_process.py ${SLURM_ARRAY_TASK_ID}

micromamba run -n GENIE  ${GENIE} \
	--genotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_${SELECTION_TYPE}_VS_${VS}_NE_${NE}.out \
	--phenotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_${SELECTION_TYPE}_VS_${VS}_NE_${NE}_phenotypes.GENIE.txt \
	--annot ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_${SELECTION_TYPE}_VS_${VS}_NE_${NE}_annotations_age_bins.txt \
	--output ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_${SELECTION_TYPE}_VS_${VS}_NE_${NE}_out_GENIE \
	--model G \
	--verbose 1 

# ### try pruning 

# micromamba run -n plink2 plink2 --bfile ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out \
#     --indep-pairwise 50 5 0.2 \
#     --out ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out

# micromamba run -n plink2 plink2 --bfile ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out \
#     --extract ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.prune.in \
#     --make-bed --out ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.pruned

#  paste \
#   <(cut -f2 ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.bim) ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.txt | \
#    grep -f ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.prune.in | cut -f2- > ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.pruned.txt


# micromamba run -n GENIE ${GENIE} \
# 	--genotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.pruned \
# 	--phenotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_phenotypes.GENIE.txt \
# 	--annot ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.pruned.txt \
# 	--output ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_neutral_out_GENIE.pruned \
# 	--model G \
# 	--verbose 1 
