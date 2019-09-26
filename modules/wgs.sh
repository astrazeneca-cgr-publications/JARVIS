#!/bin/bash 
#SBATCH -J wgs  # Set name for current job 
#SBATCH -o out.wgs  # Set job output log 
#SBATCH -e err.wgs  # Set job error log 
#SBATCH --cpus-per-task=5         # Request 5 CPUs (cores) on a single node 
#SBATCH --mem=40G          # Request amount of memory 
#SBATCH -t 48:00:0            # Request 24 hours runtime

module load libpng/1.6.23-foss-2017a


config_file=$1 
input_classes=$2 #input_classes.txt 


echo "Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.)"
# (Most time-consuming part becasue of feature extraction for each window -- GC content, mut_rate, etc.)
#./parse_all_chromosomes.sh $config_file;


echo "Perform logistic regression (common ~ all variants) to get gwRVIS scores"
#python run_full_regression.py $config_file;


echo "Convert window indexes (0-based) to real genomic coordinates"
#python convert_window_indexes_to_genomic_coords.py $config_file;




echo "Get gwRVIS distribution by genomic class"
./run_gwrvis_extraction_by_genomic_class.sh $config_file $input_classes;


echo "Compile full feature table (gwrvis, primary sequence features and regulatory features)"
python compile_full_win_feature_table.py $config_file

echo "Merge BED files by genomic class across all chromosomes"
./annotate_feat_table_w_mut_exl_genomic_class.sh $config_file $input_classes


echo "Aggregate gwRVIS scores from all chromosomes"
python aggregate_gwrvis_scores.py $config_file $input_classes;


# [Ad-hoc] post-processing: enhancers
#python process_enhancers_bed_rvis_contents.py $config_file $input_classes;

echo "Get gwRVIS distribution by genomic class across the entire genome"
python get_whole_genome_rvis_distr.py $config_file $input_classes;


#python make_ggridges_plots.py -c $config_file; # [deprecated/redundant]






# Under "scores_benchmarking/"
printf "\n\n==== Benchmarking for gwRVIS itself and against other scores (scores_benchmarking/) ===="

echo "> Run Logistic regression for gwRVIS tolerance predictive power"
python scores_benchmarking/get_gwrvis_tolerance_predictive_power.py $config_file 0
python scores_benchmarking/get_gwrvis_tolerance_predictive_power.py $config_file 1 # filtering-out gwRVIS > 0, i.e. positive selection windows



echo "> Run benchmark against clinvar/hgmd (pathogenic/benign)"
python scores_benchmarking/run_clinvar_benchmarking.py $config_file

echo "> Run benchmark against denovo-db phenotypes (cases/controls)"
#python benchmark_against_denovo_db.py $config_file

python scores_benchmarking/benchmark_vs_original_orion.py $config_file






# Under "jarvis_classification/"
printf "\n\n==== Classification with JARVIS, integrating gwRVIS and external annotations (jarvis_classification/) ===="
filter_ccds_overlapping_variants=0 # set to 1 to remove non-coding variants falling into windows that also contain coding variants
python jarvis/variant_classification/run_variant_classification.py $config_file $filter_ccds_overlapping_variants
