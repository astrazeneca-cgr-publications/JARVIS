#!/bin/bash 
#SBATCH -J wgs  # Set name for current job 
#SBATCH -o out.wgs  # Set job output log 
#SBATCH -e err.wgs  # Set job error log 
#SBATCH --cpus-per-task=5         # Request 5 CPUs (cores) on a single node 
#SBATCH --mem=40G          # Request amount of memory 
#SBATCH -t 24:00:0            # Request 24 hours runtime

module load libpng/1.6.23-foss-2017a


config_file=$1 
input_classes=$2 #input_classes.txt 

# ==================== gwRVIS core ====================

echo "Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.)"
# (Most time-consuming part becasue of feature extraction for each window -- GC content, mut_rate, etc.)
#./gwrvis_core/parse_all_chromosomes.sh $config_file;



echo "Perform logistic regression (common ~ all variants) to get gwRVIS scores"
#python gwrvis_core/run_full_regression.py $config_file;


echo "Convert window indexes (0-based) to real genomic coordinates"
#python gwrvis_core/convert_window_indexes_to_genomic_coords.py $config_file;




echo "Get gwRVIS distribution by genomic class"
#./gwrvis_core/run_gwrvis_extraction_by_genomic_class.sh $config_file $input_classes;


echo "Compile full feature table (gwrvis, primary sequence features and regulatory features)"
#python gwrvis_core/compile_full_win_feature_table.py $config_file
## <<======  Up to this point output common with BASE_OUTDIR 




echo "Merge BED files by genomic class across all chromosomes"
./gwrvis_core/annotate_feat_table_w_mut_exl_genomic_class.sh $config_file $input_classes







# ----------------------  DON'T RUN DURING DEBUG/DEV -------------------------

echo "Aggregate gwRVIS scores from all chromosomes"
#python gwrvis_core/aggregate_gwrvis_scores.py $config_file $input_classes;



# [Ad-hoc] post-processing: enhancers
##python gwrvis_core/process_enhancers_bed_rvis_contents.py $config_file $input_classes;

echo "Get gwRVIS distribution by genomic class across the entire genome"
#python gwrvis_core/get_whole_genome_rvis_distr.py $config_file $input_classes;

#python gwrvis_core/make_ggridges_plots.py -c $config_file; # [slightly redundant]






# Under "scores_benchmarking/"
printf "\n\n==== Benchmarking for gwRVIS itself and against other scores (scores_benchmarking/) ===="

echo "> Run Logistic regression for gwRVIS tolerance predictive power"
#python scores_benchmarking/get_gwrvis_tolerance_predictive_power.py $config_file 0
#python scores_benchmarking/get_gwrvis_tolerance_predictive_power.py $config_file 1 # filtering-out gwRVIS > 0, i.e. positive selection windows
#exit


echo "> Run benchmark against clinvar/hgmd (pathogenic/benign)"
#python scores_benchmarking/run_clinvar_benchmarking.py $config_file
#exit

# [To become depecreated]
echo "> Run benchmark against denovo-db phenotypes (cases/controls)"
##python scores_benchmarking/benchmark_against_denovo_db.py $config_file
#python scores_benchmarking/benchmark_vs_original_orion.py $config_file

# ------------------------------------------------------------------------------





# Update the rest of features with single-nt resolution data!
python -u jarvis/deep_learn_raw_seq/prepare_data.py $config_file



exit


# Under "jarvis_classification/"
printf "\n\n==== Classification with JARVIS, integrating gwRVIS and external annotations (jarvis/variant_classification/) ===="
filter_ccds_overlapping_variants=0 # set to 1 to remove non-coding variants falling into windows that also contain coding variants
model_type="RF"
#python jarvis/variant_classification/run_variant_classification.py $config_file $filter_ccds_overlapping_variants $model_type
exit




# Train JARVIS with structured data, sequences or both
cv_repeats=1
./submit_all_jarvis_jobs.sh $config_file $cv_repeats
exit



# ================= POST-PROCESSING ==================
# ... After all "submit_all_jarvis_jobs.sh" jobs have been complete
# # Add a monitoring script for some files that are expected to be present

# Aggregate and plot performance metrics across different genomic classes
python jarvis/performance_estimation/process_performance_metrics.py $config_file 1   #1: DL_MODELs present
