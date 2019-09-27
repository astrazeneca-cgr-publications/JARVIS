#!/bin/bash 
#SBATCH -J train_model
#SBATCH -o out.train_model
#SBATCH -e err.train_model
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16   # Request 2 CPUs (cores) on a single node
#SBATCH --mem=24G	
#SBATCH -t 24:0:0           # Request 24 hours runtime

config_file=$1
top_ratio=$2

python train_nn_model.py $config_file $top_ratio
