#!/bin/bash 
#SBATCH -J wgs  # Set name for current job 
#SBATCH -o out.wgs  # Set job output log 
#SBATCH -e err.wgs  # Set job error log 
#SBATCH --cpus-per-task=5         # Request 5 CPUs (cores) on a single node 
#SBATCH --mem=40G          # Request amount of memory 
#SBATCH -t 24:00:0            # Request 24 hours runtime

module load libpng/1.6.23-foss-2017a


config_log=$1 
input_classes=$2 #input_classes_new.txt #input_classes.txt

out_dir=`python custom_utils.py config.yaml`

echo "Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.)"
# (Most time-consuming part becasue of feature extraction for each window -- GC content, mut_rate, etc.)
#./parse_all_chromosomes.sh $config_log;


echo "Perform logistic regression (common ~ all variants) to get gwRVIS scores"
#python run_full_regression.py $config_log;


echo "Convert window indexes (0-based) to real genomic coordinates"
#python convert_window_indexes_to_genomic_coords.py $config_log;





echo "Get gwRVIS distribution by genomic class"
#./run_gwrvis_extraction_by_genomic_class.sh $config_log $input_classes;


echo "Compile full feature table (gwrvis, primary sequence features and regulatory features)"
#python compile_full_win_feature_table.py $config_log

echo "Merge BED files by genomic class across all chromosomes"
./annotate_feat_table_w_mut_exl_genomic_class.sh $out_dir $input_classes
exit

echo "Aggregate gwRVIS scores from all chromosomes"
python aggregate_gwrvis_scores.py $config_log $input_classes;


# [Deprecated] post-processing: enhancers
#python process_enhancers_bed_rvis_contents.py $config_log $input_classes;

echo "Get gwRVIS distribution by genomic class across the entire genome"
python get_whole_genome_rvis_distr.py $config_log $input_classes;


#python make_ggridges_plots.py -c $config_log; 


echo "Run Logistic regression for gwRVIS predictive power"
python ml_modules/run_gwrvis_logistic_regression.py $config_log 0
python ml_modules/run_gwrvis_logistic_regression.py $config_log 1 # filtering-out gwRVIS > 0, i.e. positive selection windows


echo "Benchmarking of gwRVIS against other scores"
cd scores_benchmarking

echo "Run benchmark against clinvar (pathogenic/benign)"
python run_clinvar_benchmarking.py ../$config_log

echo "Run benchmakr against denovo-db phenotypes (cases/controls)"
#python benchmark_against_denovo_db.py ../$config_log

# DEBUG: relative and absolute paths
#python benchmark_gwrvis_vs_original_orion.py ../$config_log
