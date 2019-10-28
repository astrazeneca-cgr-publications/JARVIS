#!/bin/bash

config_file=$1
get_new_fixed_cv_batch=0 # 1: create new fixed CV batch, 0: otherwise



#declare -a genomic_classes=("intergenic" "utr" "intergenic,utr" "lincrna" "intergenic,utr,lincrna,ucne,vista" "intron" "ccds")
declare -a genomic_classes=("intergenic")

# > Train for "structured" features first and get a fixed CV batch set for use with each genomic class
if [ $get_new_fixed_cv_batch == 1 ]; then
	echo "Training JARVIS with structured features..."

	use_fixed_cv_batches=0
	in_features="structured"
	for gen_classes in "${genomic_classes[@]}"; do

		job_name="${in_features}_${gen_classes}"
		echo "Job: $job_name"

		if [ $in_features == 'sequences' ] || [ $in_features == 'both' ]; then
			ncores=8
		else
			ncores=4
		fi

		# ClinVar pathogenic set
		sbatch -o "logs/clinvar.${job_name}.out" -e "logs/clinvar.${job_name}.err" --mem-per-cpu=24G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $in_features $gen_classes $use_fixed_cv_batches

		# HGMD pathogenic set
		#sbatch -o "logs/hgmd.${job_name}.out" -e "logs/hgmd.${job_name}.err" --mem-per-cpu=36G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $in_features $gen_classes $use_fixed_cv_batches
	done
fi




# > Train for "sequences" and "both" features using fixed CV batches built earlier while training for "structured" features 
declare -a input_features=("sequences" "both")
use_fixed_cv_batches=1

for in_features in "${input_features[@]}"; do
	echo "Training JARVIS with $in_features features..."

	for gen_classes in "${genomic_classes[@]}"; do

		job_name="${in_features}_${gen_classes}"
		echo "Job: $job_name"

		if [ $in_features == 'sequences' ] || [ $in_features == 'both' ]; then
			ncores=8
		else
			ncores=4
		fi

		# ClinVar pathogenic set
		sbatch -o "logs/clinvar.${job_name}.out" -e "logs/clinvar.${job_name}.err" --mem-per-cpu=24G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $in_features $gen_classes $use_fixed_cv_batches

		# HGMD pathogenic set
		#sbatch -o "logs/hgmd.${job_name}.out" -e "logs/hgmd.${job_name}.err" --mem-per-cpu=36G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $in_features $gen_classes $use_fixed_cv_batches

	done
done
