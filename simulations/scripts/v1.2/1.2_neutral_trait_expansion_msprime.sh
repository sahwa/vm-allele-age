#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J v1.2.1_neutral_trait_expansion_msprime
#SBATCH -o v1.2.1_neutral_trait_expansion_msprime.%A.out
#SBATCH -e v1.2.1_neutral_trait_expansion_msprime.%A.err
#SBATCH -p short
#SBATCH -c 4

micromamba run -n slim python v1.2.1_neutral_trait_expansion_msprime.py
