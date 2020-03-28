#!/bin/bash

config_file=$1


if [ "$#" -ne 1 ]; then
	echo -e "[Error] Insufficient args. Re-run with 2 input args: [config_file]\n"
	exit
fi

declare -a genomic_classes=("intergenic,utr,lincrna,ucne,vista") 
cv_repeats=1


for gen_classes in "${genomic_classes[@]}"; do

	job_name="${gen_classes}.jarvis"
	echo -e "\nJob: $job_name"

	ncores=1
	mem="40G"
	t="24:0:0"

	sbatch -o "logs/pred.${job_name}.out" -e "logs/pred.${job_name}.err" --time=$t --mem-per-cpu=${mem} --cpus-per-task=${ncores} ./jarvis/deep_learn_raw_seq/predict_on_test.sh $config_file $gen_classes $cv_repeats
done
