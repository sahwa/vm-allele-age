#!/bin/bash
#SBATCH -A visscher-wray.prj
#SBATCH -J 1.1_neutral_trait.v2
#SBATCH -o 1.1_neutral_trait.v2.%A_%a.out
#SBATCH -e 1.1_neutral_trait.v2.%A_%a.err
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 12:00:00

slim -d END_TICK=100000 -d outfile="'1.1_neutral_trait.n_100000.v2'" 1.1_neutral_trait.n_100000.slim 

