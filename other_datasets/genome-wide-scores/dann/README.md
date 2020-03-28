# Data preparation
sbatch ./preprocess_dann_raw_scores.sh



# Get DANN subset BED files for different types of variant sets (either pathogenic or benign, from various resources)
sbatch ./prepare_variant_specific_files.sh
