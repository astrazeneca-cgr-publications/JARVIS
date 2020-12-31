#!/bin/bash 
#SBATCH -t 48:0:0           
##SBATCH --mem=150G
#SBATCH --mem=200G  

## CPU
##SBATCH --cpus-per-task=20

## GPU
#SBATCH --partition=gpu
#SBATCH --gres=gpu:volta:4



config_file=$1
input_features=$2
chrom=$3

module load libpng/1.6.23-foss-2017a


python -u jarvis/deep_learn_raw_seq/predict_genome_wide_jarvis.py $config_file $input_features $chrom
