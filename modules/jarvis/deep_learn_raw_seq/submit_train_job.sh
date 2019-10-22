#!/bin/bash 
#SBATCH -t 24:0:0           # Request 24 hours runtime
#SBATCH --mem-per-cpu=24G	# default values - may be overriden during sbatch call
#SBATCH --cpus-per-task=4	# default values - may be overriden during sbatch call

config_file=$1
input_features=$2	# structured, sequences or both
genomic_classes=$3  	# comma-separated if more than 1
use_fixed_cv_batches=$4


python jarvis/deep_learn_raw_seq/train_nn_model.py $config_file $input_features $genomic_classes $use_fixed_cv_batches
