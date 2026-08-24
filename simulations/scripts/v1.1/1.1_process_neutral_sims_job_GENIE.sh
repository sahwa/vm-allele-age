#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_process
#SBATCH -o 1.1_process.%A.out
#SBATCH -e 1.1_process.%A.err
#SBATCH -p short
#SBATCH -c 5
#SBATCH --mem=48G
#SBATCH -t 08:00:00

VERSION=1.1

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}

VCF=${DATA}/${VERSION}_neutral_out.vcf

plink2 --vcf ${VCF} \
	--make-bed \
	--max-alleles 2 \
	--min-alleles 2 \
	--out ${DATA}/${VERSION}_neutral_out


GENIE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/GENIE/build/GENIE

${GENIE} \
	--genotype ${DATA}/${VERSION}_neutral_out \
	--phenotype ${DATA}/1.1_phenotypes.GENIE.txt \
	--annot ${DATA}/1.1_annotations_age_bins.txt \
	--output ${DATA}/1.1_neutral_out_GENIE \
	--model G

### testing using just one bin ###
awk '{print 1}' ${DATA}/${VERSION}_neutral_out.bim > ${DATA}/single_annot.txt
${GENIE} --genotype ${DATA}/${VERSION}_neutral_out \
         --phenotype ${DATA}/1.1_phenotypes.GENIE.txt \
         --annot ${DATA}/single_annot.txt \
         --output ${DATA}/1.1_test_single \
         --model G
		 