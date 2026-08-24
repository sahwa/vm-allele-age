#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_neutral_trait.submit_slim
#SBATCH -o 1.1_neutral_trait.submit_slim.%A_%a.out
#SBATCH -e 1.1_neutral_trait.submit_slim.%A_%a.err
#SBATCH -p short
#SBATCH -c 6

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.1

cd ${DATA}

bins="$(ls 1.1_bin_*.bed | sed 's/.bed//g')"

for bin in `ls 1.1_bin_*.bed | sed 's/.bed//g'`; do
    # plink2 --vcf ${DATA}/1.0_neutral_out.vcf.gz --extract ${bin} --make-bed --out ${DATA}/1.0_neutral_out.${bin_name}
    gcta --bfile ${DATA}/${bin} --make-grm-bin --out ${DATA}/${bin}
    gcta64 --grm ${DATA}/${bin} --grm-summary --out ${DATA}/${bin}
done


ls 1.1_bin_*_*.grm.bin | \
    cut -d'.' -f1-2 | \
    sort -V |
    > 1.1_neutral_out.variant_info_age.grm.bin.list

gcta64 --reml --mgrm 1.1_neutral_out.variant_info_age.grm.bin.list --pheno ${DATA}/1.1_phenotypes.txt --out 1.1_neutral_out.phenotypes
