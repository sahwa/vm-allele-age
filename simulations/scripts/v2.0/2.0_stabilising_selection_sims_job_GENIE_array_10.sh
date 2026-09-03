#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 2.0_process_array
#SBATCH -o logs/2.0_process.%A_%a.out
#SBATCH -e logs/2.0_process.%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH -a 1-3

VERSION=2.0
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v${VERSION}
SELECTION_TYPE="stabilising_selection"
VS=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${SCRIPTS}/2.0_stabilising_selection_VS_sweep.txt)
# VS=5
NE=10000
PARAMS=VS_${VS}_NE_${NE}
FILE_STEM=${VERSION}_${SELECTION_TYPE}_${PARAMS}
GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}

for REP in {1..10}; do

 DATA_REP=${DATA}/replicates/${PARAMS}/${REP}

 micromamba run -n slim python ${SCRIPTS}/2.0_stabilising_selection_process.py \
  --tree ${DATA}/2.0_stabilising_selection_VS_${VS}_NE_${NE}.trees \
  --vs ${VS} \
  --ne ${NE} \
  --rep ${REP}

 ${GENIE} \
  --genotype ${DATA_REP}/${FILE_STEM}.out \
  --phenotype ${DATA_REP}/${FILE_STEM}_phenotypes.GENIE.txt \
  --annot ${DATA_REP}/${FILE_STEM}_annotations_age_bins.txt \
  --output ${DATA_REP}/${FILE_STEM}_out_GENIE \
  --model G \
  --verbose 1 \
  --num-jack 20
done

# ---------------------------------------------------------------
# LD pruning (commented out until decided whether it's needed for
# the selection sims - the neutral v1.1 case showed it relocates
# rather than removes bias, so worth checking here separately
# before turning this back on)
# ---------------------------------------------------------------

# micromamba run -n plink2 plink2 --bfile ${DATA_REP}/${FILE_STEM}.out \
#     --indep-pairwise 50 5 0.2 \
#     --out ${DATA_REP}/${FILE_STEM}.out
#
# micromamba run -n plink2 plink2 --bfile ${DATA_REP}/${FILE_STEM}.out \
#     --extract ${DATA_REP}/${FILE_STEM}.out.prune.in \
#     --make-bed --out ${DATA_REP}/${FILE_STEM}.out.pruned
#
# awk 'NR==FNR {keep[$1]; next} ($1 in keep)' \
#     ${DATA_REP}/${FILE_STEM}.out.prune.in \
#     <(paste <(cut -f2 ${DATA_REP}/${FILE_STEM}.out.bim) \
#             ${DATA_REP}/${FILE_STEM}_annotations_age_bins.txt) \
#     | cut -f2- > ${DATA_REP}/${FILE_STEM}_annotations_age_bins.pruned.txt
#
# test $(wc -l < ${DATA_REP}/${FILE_STEM}_annotations_age_bins.pruned.txt) -eq \
#      $(wc -l < ${DATA_REP}/${FILE_STEM}.out.pruned.bim) \
#      || { echo "ANNOTATION/BIM MISMATCH rep${REP}"; exit 1; }
#
# micromamba run -n GENIE ${GENIE} \
#     --genotype ${DATA_REP}/${FILE_STEM}.out.pruned \
#     --phenotype ${DATA_REP}/${FILE_STEM}_phenotypes.GENIE.txt \
#     --annot ${DATA_REP}/${FILE_STEM}_annotations_age_bins.pruned.txt \
#     --output ${DATA_REP}/${FILE_STEM}_out_GENIE.pruned \
#     --model G \
#     --verbose 1