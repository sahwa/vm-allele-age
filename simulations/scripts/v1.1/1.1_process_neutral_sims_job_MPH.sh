#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_process_array
#SBATCH -o 1.1_process.%A_%a.out
#SBATCH -e 1.1_process.%A_%a.err
#SBATCH -p short
#SBATCH -c 6
#SBATCH -t 04:00:00
#SBATCH -a 1-10

REP=${SLURM_ARRAY_TASK_ID}  # <- ADD THIS

VERSION=1.1

DATA=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v${VERSION}
SCRIPTS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/scripts/v1.1
PROGRAMS=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/programs/

GENIE=${PROGRAMS}/GENIE/build/GENIE
MPH=${PROGRAMS}/mph/mph

# micromamba run -n slim python ${SCRIPTS}/1.1_process_neutral_sims.py ${REP}

for bin in $(ls ${DATA}/rep${REP}/1.1_bin_*-*.bed | sed 's/\.bed$//'); do
    stem="$(basename ${bin})"577718

    micromamba run -n MPH ${MPH} --bfile ${bin} --make_grm \
        --snp_info_file ${DATA}/rep${REP}/1.1_snp_info_mph.csv \
        --snp_weight_name w_unstd \
        --snp_genotype_coding c0,c1,c2 \
        --output_file ${DATA}/rep${REP}/${stem}_MPH_unstd \
        --num_threads 6

    micromamba run -n MPH ${MPH} --bfile ${bin} --make_grm \
        --snp_info_file ${DATA}/rep${REP}/1.1_snp_info_mph.csv \
        --snp_weight_name w_std \
        --snp_genotype_coding c0,c1,c2 \
        --output_file ${DATA}/rep${REP}/${stem}_MPH_std \
        --num_threads 6
done


# #### then make the phenotype file ####
awk -F',' 'NR>1 {print $2","$3}' ${DATA}/rep${REP}/1.1_phenotypes.csv \
    | sed '1i id,PHENO' > ${DATA}/rep${REP}/1.1_MPH_pheno.txt

#### run multi-component REML for standardised GRMs ####

ls ${DATA}/rep${REP}/*_MPH_std*.grm.bin | sed 's/\.grm\.bin$//' > ${DATA}/rep${REP}/1.1_MPH_std.grm.list
micromamba run -n MPH ${MPH} \
    --grm_list ${DATA}/rep${REP}/1.1_MPH_std.grm.list \
    --phenotype_file ${DATA}/rep${REP}/1.1_MPH_pheno.txt \
    --trait PHENO \
    --reml \
    --heritability 0.5 \
    --seed 42 \
    --output_file ${DATA}/rep${REP}/1.1_MPH_std \
    --num_threads 6 

#### run multi-component REML for unstandardised GRMs ####
ls ${DATA}/rep${REP}/*_MPH_unstd*.grm.bin | sed 's/\.grm\.bin$//' > ${DATA}/rep${REP}/1.1_MPH_unstd.grm.list
micromamba run -n MPH ${MPH} \
    --grm_list ${DATA}/rep${REP}/1.1_MPH_unstd.grm.list \
    --phenotype_file ${DATA}/rep${REP}/1.1_MPH_pheno.txt \
    --trait PHENO \
    --reml \
    --heritability 0.5 \
    --seed 42 \
    --output_file ${DATA}/rep${REP}/1.1_MPH_unstd \
    --num_threads 6
    
