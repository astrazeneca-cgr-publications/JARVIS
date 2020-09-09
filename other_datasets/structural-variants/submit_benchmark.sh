#!/bin/bash
#SBATCH -t 24:0:0
#SBATCH -n 1
#SBATCH --mem=150G

score=$1

python -u benchmark_scores.py $score
