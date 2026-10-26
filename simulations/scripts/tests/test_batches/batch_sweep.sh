#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J batch_sweep
#SBATCH -o batch_sweep.%A.out
#SBATCH -e batch_sweep.%A.err
#SBATCH -p short
#SBATCH -c 4

TREE=/exafs1/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/v1.2/n500000_L1e+08/sim.trees
OUT=/well/visscher-wray/users/uwu199/projects/vm-allele-age/simulations/data/tests/batch_test

micromamba run -n slim python batch_sweep.py \
    --trees ${TREE} \
    --out ${OUT} \
   --n-batches 1 2 5 10 20 50 100 200


