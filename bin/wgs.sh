#!/bin/bash 
#SBATCH -J wgs  # Set name for current job 
#SBATCH -o out.wgs  # Set job output log 
#SBATCH -e err.wgs  # Set job error log 
#SBATCH --cpus-per-task=5         # Request 5 CPUs (cores) on a single node 
#SBATCH --mem=40G          # Request amount of memory 
#SBATCH -t 24:00:0            # Request 24 hours runtime

module load libpng/1.6.23-foss-2017a


config_log=$1 
input_classes=input_classes.txt

# Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.)
# (Most time-consuming part becasue of feature extraction for each window -- GC content, mut_rate, etc.)
#./parse_all_chromosomes.sh $config_log;
#exit

# Perform logistic regression (common ~ all variants) to get gwRVIS scores
python run_full_regression.py $config_log;


# Convert window indexes (0-based) to real genomic coordinates
python convert_window_indexes_to_genomic_coords.py $config_log;


# Get gwRVIS distribution by genomic class
./run_gwrvis_extraction_by_genomic_class.sh $config_log $input_classes;


# Aggregate gwRVIS scores from all chromosomes
python aggregate_gwrvis_scores.py $config_log $input_classes;


# [Deprecated] post-processing: enhancers
#python process_enhancers_bed_rvis_contents.py config.log input_classes.txt;

python get_whole_genome_rvis_distr.py $config_log $input_classes;


#python make_ggridges_plots.py -c $config_log; 


# Run Logistic regression for gwRVIS predictive power
python ml_modules/run_gwrvis_logistic_regression.py $config_log 0
python ml_modules/run_gwrvis_logistic_regression.py $config_log 1 # filtering-out gwRVIS > 0, i.e. positive selection windows


out_dir=`python custom_utils.py config.yaml`
