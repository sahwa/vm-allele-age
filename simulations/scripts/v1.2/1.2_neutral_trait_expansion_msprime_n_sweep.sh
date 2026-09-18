#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J v1.2_neutral_trait_expansion_msprime_n_sweep
#SBATCH -o v1.2_neutral_trait_expansion_msprime_n_sweep.%A.out
#SBATCH -e v1.2_neutral_trait_expansion_msprime_n_sweep.%A.err
#SBATCH -p long
#SBATCH -c 4
#SBATCH -a 1-10


line=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" 1.2_neutral_trait_expansion_msprime_n_sweep.params)
N=$(echo $line | cut -d' ' -f1)
L=$(echo $line | cut -d' ' -f2)

micromamba run -n slim python 1.2_neutral_trait_expansion_msprime_n_sweep.py ${N} ${L}
