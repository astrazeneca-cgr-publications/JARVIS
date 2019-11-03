#!/bin/bash

config_file=$1
cv_repeats=$2


if [ "$#" -ne 2 ]; then
	echo -e "[Error] Insufficient args. Re-run with 2 input args: [config_file] [cv_repeats]\n"
	exit
fi


#declare -a genomic_classes=("intergenic" "utr" "intergenic,utr" "lincrna" "intergenic,utr,lincrna,ucne,vista") 
#declare -a genomic_classes=("ccds" "intron")
declare -a genomic_classes=("intron")


for gen_classes in "${genomic_classes[@]}"; do

	job_name="${gen_classes}.jarvis"
	echo -e "\nJob: $job_name"

	ncores=4
	mem="24G"
	t="24:0:0"
	if [ "$gen_classes" == "ccds" ] || [ "$gen_classes" == "intron" ]; then
		ncores=10
		mem="24G"
		t="48:0:0"
	fi

	# ClinVar pathogenic set
	sbatch -o "logs/clinvar.${job_name}.out" -e "logs/clinvar.${job_name}.err" --time=$t --mem-per-cpu=${mem} --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $gen_classes $cv_repeats

	# HGMD pathogenic set
	#sbatch -o "logs/hgmd.${job_name}.out" -e "logs/hgmd.${job_name}.err" --time=$t --mem-per-cpu=36G --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $gen_classes $cv_repeats
done

