#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_process
#SBATCH -o 1.1_process.%A.out
#SBATCH -e 1.1_process.%A.err
#SBATCH -p short
#SBATCH -c 5
#SBATCH --mem=48G
#SBATCH -t 08:00:00

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1

python 1.1_process_neutral_sims.py

bins="$(ls ${DATA}/1.0_neutral_out.variant_info_age_*csv)"

for bin in 1.0_neutral_out.variant_info_age_*csv; do
    bin_name=$(echo $bin | cut -d'.' -f3)
    plink2 --vcf ${DATA}/1.0_neutral_out.vcf.gz --extract ${bin} --make-bed --out ${DATA}/1.0_neutral_out.${bin_name}
    gcta --bfile ${DATA}/1.0_neutral_out.${bin_name} --make-grm-bin --out ${DATA}/1.0_neutral_out.${bin_name}
done

ls 1.0_neutral_out.variant_info_age_*.grm.bin | \
cut -d'.' -f1-3 | \
sort -V > 1.0_neutral_out.variant_info_age.grm.bin.list

gcta64 --reml --mgrm ${DATA}/1.0_neutral_out.variant_info_age.grm.bin.list --pheno ${DATA}/1.0_neutral_out.phenotypes --out ${DATA}/1.0_neutral_out.phenotypes

