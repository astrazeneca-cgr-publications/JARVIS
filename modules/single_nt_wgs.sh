#!/bin/bash 
#SBATCH --cpus-per-task=3
#SBATCH --mem=35G
#SBATCH -t 24:00:0

module load libpng/1.6.23-foss-2017a

config_file=$1 
single_nt_offset=$2   

chr="-1" # deprecated: 
# cannot run for a single chromosome, need to have the common/all variant counts across all chromosomes to run the regression - can do that separately for each nt-offset though





# ==================== gwRVIS core ====================
function get_gwrvis_per_offset {

	single_nt_offset=$1

	echo "Record features across sliding genomic windows (common/all variants)"
	# (Most time-consuming part becasue of feature extraction for each window -- GC content, mut_rate, etc.)
	./gwrvis_core/single_nt_parse_all_chromosomes.sh $config_file $single_nt_offset $chr;



	#echo "Perform logistic regression (common ~ all variants) to get gwRVIS scores"
	python gwrvis_core/single_nt_run_full_regression.py $config_file $single_nt_offset $chr;



	#echo "Convert window indexes (0-based) to real genomic coordinates"
	python gwrvis_core/single_nt_convert_window_indexes_to_genomic_coords.py $config_file $single_nt_offset $chr;

}


get_gwrvis_per_offset $single_nt_offset






