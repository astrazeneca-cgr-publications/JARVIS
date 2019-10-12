#!/bin/bash


declare -a input_features=("both") #(structured sequences both)
declare -a genomic_classes=("intergenic" "utr" "intergenic,utr" "lincrna" "intergenic,utr,lincrna,ucne,vista" "intron" "ccds")



for in_features in "${input_features[@]}"; do
	for gen_classes in "${genomic_classes[@]}"; do

		job_name="${in_features}_${gen_classes}"
		echo "Job: $job_name"

		if [ $in_features == 'sequences' ] || [ $in_features == 'both' ]; then
			ncores=8
		else
			ncores=4
		fi

		# ClinVar pathogenic set
		sbatch -o "logs/clinvar.${job_name}.out" -e "logs/clinvar.${job_name}.err" --mem-per-cpu=24G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh conf/config.yaml $in_features $gen_classes

		# HGMD pathogenic set
		sbatch -o "logs/hgmd.${job_name}.out" -e "logs/hgmd.${job_name}.err" --mem-per-cpu=36G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh conf/config.hgmd.yaml $in_features $gen_classes

	done
done
