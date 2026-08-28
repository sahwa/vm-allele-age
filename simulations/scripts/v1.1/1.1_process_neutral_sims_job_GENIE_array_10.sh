#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_process_array
#SBATCH -o 1.1_process.%A_%a.out
#SBATCH -e 1.1_process.%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH -a 1-10

VERSION=1.1

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v1.1

GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

#python ${SCRIPTS}/1.1_process_neutral_sims.py ${SLURM_ARRAY_TASK_ID}

# ${GENIE} \
# 	--genotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out \
# 	--phenotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_phenotypes.GENIE.txt \
# 	--annot ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.txt \
# 	--output ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_neutral_out_GENIE \
# 	--model G \
# 	--verbose 1 

### try pruning 

micromamba run -n plink2 plink2 --bfile ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out \
    --indep-pairwise 50 5 0.2 \
    --out ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out

micromamba run -n plink2 plink2 --bfile ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out \
    --extract ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.prune.in \
    --make-bed --out ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.pruned

 paste \
  <(cut -f2 ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.bim) ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.txt | \
   grep -f ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.prune.in | cut -f2- > ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.pruned.txt


micromamba run -n GENIE ${GENIE} \
	--genotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/${VERSION}_neutral_out.pruned \
	--phenotype ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_phenotypes.GENIE.txt \
	--annot ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_annotations_age_bins.pruned.txt \
	--output ${DATA}/rep${SLURM_ARRAY_TASK_ID}/1.1_neutral_out_GENIE.pruned \
	--model G \
	--verbose 1 
