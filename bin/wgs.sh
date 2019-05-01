#!/bin/bash 
#SBATCH -J wgs  # Set name for current job 
#SBATCH -o out.wgs  # Set job output log 
#SBATCH -e err.wgs  # Set job error log 
#SBATCH --cpus-per-task=5         # Request 23 CPUs (cores) on a single node 
#SBATCH --mem=40000          # Request amount of memory 
#SBATCH -t 0:30:0            # Request 24 hours runtime

config_log=config.log;
input_classes=input_classes.txt;

# Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.)
./parse_all_chromosomes.sh $config_log;


# Perform logistic regression (common ~ all variants) to get gwRVIS scores
python run_full_regression.py $config_log;


# Convert window indexes (0-based) to real genomic coordinates
python convert_window_indexes_to_genomic_coords.py $config_log;


# Get gwRVIS distribution by genomic class
./run_gwrvis_extraction_by_genomic_class.sh $config_log $input_classes;


# Aggregate gwRVIS scores from all chromosomes
python aggregate_gwrvis_scores.py $config_log $input_classes;

# ---------------


# post-processing: enhancers
#python process_enhancers_bed_rvis_contents.py config.log input_classes.txt;


# Available option: filter out rvis scores from intergenic regions with no variation
python whole_genome_rvis_distr.py $config_log $input_classes; # 0;


python make_whole_genome_boxplots.py $config_log;


python post_process_results.py -c $config_log; 


