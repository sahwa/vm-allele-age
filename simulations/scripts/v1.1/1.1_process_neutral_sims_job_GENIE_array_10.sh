#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_process_array
#SBATCH -o 1.1_process.%A_%a.out
#SBATCH -e 1.1_process.%A_%a.err
#SBATCH -p short
#SBATCH -c 3
#SBATCH -a 1-100

VERSION=1.1

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v1.1

GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

REP_BASE=${DATA}/replicates/rep${SLURM_ARRAY_TASK_ID}


micromamba run -n slim python ${SCRIPTS}/1.1_process_neutral_sims.py ${SLURM_ARRAY_TASK_ID}

micromamba run -n GENIE  ${GENIE} \
 --genotype ${REP_BASE}/${VERSION}_neutral_out \
 --phenotype ${REP_BASE}/${VERSION}_phenotypes.GENIE.txt \
 --annot ${REP_BASE}/1.1_annotations_age_bins.txt \
 --output ${REP_BASE}/1.1_neutral_out_GENIE \
 --model G \
 --verbose 1 \
 --nthreads 2

### try pruning 

micromamba run -n plink2 plink2 --bfile ${REP_BASE}/${VERSION}_neutral_out \
    --indep-pairwise 50 5 0.2 \
    --out ${REP_BASE}/${VERSION}_neutral_out

micromamba run -n plink2 plink2 --bfile ${REP_BASE}/${VERSION}_neutral_out \
    --extract ${REP_BASE}/${VERSION}_neutral_out.prune.in \
    --make-bed --out ${REP_BASE}/${VERSION}_neutral_out.pruned

awk 'NR==FNR {keep[$1]; next} ($1 in keep)' \
    ${REP_BASE}/${VERSION}_neutral_out.prune.in \
    <(paste <(cut -f2 ${REP_BASE}/${VERSION}_neutral_out.bim) \
            ${REP_BASE}/1.1_annotations_age_bins.txt) \
    | cut -f2- > ${REP_BASE}/1.1_annotations_age_bins.pruned.txt


micromamba run -n GENIE ${GENIE} \
	--genotype ${REP_BASE}/${VERSION}_neutral_out.pruned \
	--phenotype ${REP_BASE}/1.1_phenotypes.GENIE.txt \
	--annot ${REP_BASE}/1.1_annotations_age_bins.pruned.txt \
	--output ${REP_BASE}/1.1_neutral_out_GENIE.pruned \
	--model G \
	--verbose 1 \
	--nthreads 2
