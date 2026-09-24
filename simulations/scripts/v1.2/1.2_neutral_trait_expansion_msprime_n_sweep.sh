#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J v1.2_neutral_trait_expansion_msprime_n_sweep
#SBATCH -o v1.2_neutral_trait_expansion_msprime_n_sweep.%A_%a.out
#SBATCH -e v1.2_neutral_trait_expansion_msprime_n_sweep.%A_%a.err
#SBATCH -p long
#SBATCH -c 4
#SBATCH -a 1-10


BASE=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2
REF_JSON=${BASE}/n500000_L1e+08/reference_constants.json

line=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" 1.2_neutral_trait_expansion_msprime_n_sweep.params)
N=$(echo $line | cut -d' ' -f1)
L=$(echo $line | cut -d' ' -f2)

### this actually runs the sims ##
# micromamba run -n slim python 1.2_neutral_trait_expansion_msprime_n_sweep.py ${N} ${L}

## this actually sweeps over the sim output and i) sweeps over the S params and then ii)
## writes the outpt

N="$(echo $N | sed 's/_//g')"
L="$(printf "%.0e\n" $L)"

micromamba run -n slim python 1.2_neutral_trait_expansion_msprime_n_sweep_phenotypes.py \
	--tree ${BASE}/n${N}_L${L}/sim.trees \
	--ref-json ${REF_JSON} \
	--out ${BASE}/n${N}_L${L}
