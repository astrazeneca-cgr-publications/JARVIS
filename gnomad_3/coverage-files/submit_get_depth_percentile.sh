#!/bin/bash
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:0:0

python get_depth_percentile.py
