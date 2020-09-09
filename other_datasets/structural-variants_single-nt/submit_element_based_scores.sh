#!/bin/bash
#SBATCH --mem=16G
#SBATCH -n 4
#SBATCH -t 24:00:00

score=$1
class=$2

python -u sv_benchmark_element_based.py $score $class
