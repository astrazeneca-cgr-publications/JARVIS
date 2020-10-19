#!/bin/bash

config_file=$1
cv_repeats=$2


if [ "$#" -ne 2 ]; then
	echo -e "[Error] Insufficient args. Re-run with 2 input args: [config_file] [cv_repeats]\n"
	exit
fi


#declare -a genomic_classes=("intergenic" "utr" "lincrna" "intergenic,utr,lincrna,ucne,vista") 
declare -a genomic_classes=("intergenic,utr,lincrna,ucne,vista") 
#declare -a genomic_classes=("intergenic")
#declare -a genomic_classes=("ccds" "intron")



# [REDUNDANT] 
# -- now part of wgs.sh - to be called before RF/GB classifications
# -- Create jarvis_data.pkl to be used in all other runs without conflicts
#python -u jarvis/deep_learn_raw_seq/prepare_data.py $config_file



for gen_classes in "${genomic_classes[@]}"; do

	job_name="${gen_classes}.jarvis"
	echo -e "\nJob: $job_name"

	ncores=4
	mem="24G"
	t="24:0:0"
	if [ "$gen_classes" == "ccds" ] || [ "$gen_classes" == "intron" ]; then
		ncores=4
		mem="24G"
		t="48:0:0"
	fi

	# ClinVar pathogenic set
	sbatch -o "logs/clinvar.${job_name}.out" -e "logs/clinvar.${job_name}.err" --time=$t --mem-per-cpu=${mem} --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/submit_train_job.sh $config_file $gen_classes $cv_repeats

done

