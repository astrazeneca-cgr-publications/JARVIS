#!/bin/bash
#SBATCH -o out.gwrvis.plot_distr
#SBATCH --time=24:00:00
#SBATCH --mem=100G

score=$1
genomic_class=$2

python -u plot_gen_pop_vs_sv_distr.py $score $genomic_class
