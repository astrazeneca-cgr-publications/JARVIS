#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G 
#SBATCH --time=24:0:0

python -u prepare_wins_by_intol_and_conserv.py
