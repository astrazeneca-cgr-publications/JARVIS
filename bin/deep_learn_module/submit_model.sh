#!/bin/bash 
#SBATCH -J rnn_train_model
#SBATCH -o out.rnn_train_model
#SBATCH -e err.rnn_train_model
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8   # Request 2 CPUs (cores) on a single node
#SBATCH --mem=48G	
#SBATCH -t 24:0:0           # Request 24 hours runtime

top_ratio=$1

python train_nn_model.py ../config.yaml $top_ratio
