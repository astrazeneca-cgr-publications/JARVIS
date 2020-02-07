#!/bin/bash 
#SBATCH -t 24:0:0           # Request 24 hours runtime

## CPU
#SBATCH --mem-per-cpu=50G	# default values - may be overriden during sbatch call
#SBATCH --cpus-per-task=1	# default values - may be overriden during sbatch call

config_file=$1
genomic_classes=$2  	# comma-separated if more than 1
cv_repeats=1

module load libpng/1.6.23-foss-2017a


# ----------- Structured -----------
#use_fixed_cv_batches=0
#echo -e "\nTraining JARVIS for $genomic_classes using 'structured' features..."
#python -u jarvis/deep_learn_raw_seq/train_nn_model.py $config_file structured $genomic_classes $use_fixed_cv_batches $cv_repeats


# -----------  Sequences -----------
#use_fixed_cv_batches=1
#echo -e "\nTraining JARVIS for $genomic_classes using 'sequences' features..."
#python -u jarvis/deep_learn_raw_seq/train_nn_model.py $config_file sequences $genomic_classes $use_fixed_cv_batches $cv_repeats


# -----------  Both -----------
use_fixed_cv_batches=0
echo -e "\nTraining JARVIS for $genomic_classes using 'both' features..."
python -u jarvis/deep_learn_raw_seq/train_nn_model.py $config_file both $genomic_classes $use_fixed_cv_batches $cv_repeats
