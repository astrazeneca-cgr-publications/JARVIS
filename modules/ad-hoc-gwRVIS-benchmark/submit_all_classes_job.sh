#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -o boxplots.linc_intergenic.out
#SBATCH -n 1
#SBATCH --mem=200G    #400G

win_len=3000
python plot_gwrvis_distr_by_genomic_class.py $win_len all
