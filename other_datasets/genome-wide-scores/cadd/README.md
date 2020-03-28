# Data preparation
sbatch ./preprocess_cadd_raw_scores.sh

# Get CADD subset BED files for different types of variant sets (either pathogenic or benign, from various resources)
sbatch ./prepare_variant_specific_files.sh
