#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -o W3000.logit_regression.out
##SBATCH -o W3000.logit_regression.discard_positives.out
#SBATCH -n 1
#SBATCH --mem=300G

win_len=$1     # default: 3000
discard_positives=$2	# default: 0
python -u get_gwrvis_tolerance_predictive_power.py $win_len $discard_positives
