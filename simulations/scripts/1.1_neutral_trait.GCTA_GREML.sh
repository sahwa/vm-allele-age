vim#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_neutral_trait.submit_slim
#SBATCH -o 1.1_neutral_trait.submit_slim.%A_%a.out
#SBATCH -e 1.1_neutral_trait.submit_slim.%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=32G

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data

bins="$(ls ${DATA}/1.0_neutral_out.variant_info_age_*csv)"

for bin in 1.0_neutral_out.variant_info_age_*csv; do
    bin_name=$(echo $bin | cut -d'.' -f3)
    plink2 --vcf ${DATA}/1.0_neutral_out.vcf.gz --extract ${bin} --make-bed --out ${DATA}/1.0_neutral_out.${bin_name}
    gcta --bfile ${DATA}/1.0_neutral_out.${bin_name} --make-grm-bin --out ${DATA}/1.0_neutral_out.${bin_name}
done

ls 1.0_neutral_out.variant_info_age_*.grm.bin | \
cut -d'.' -f1-3 | \
sort -V > 1.0_neutral_out.variant_info_age.grm.bin.list

gcta64 --reml --mgrm 1.0_neutral_out.variant_info_age.grm.bin.list --pheno 1.0_neutral_out.phenotypes --out 1.0_neutral_out
