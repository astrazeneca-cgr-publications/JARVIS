#!/bin/bash

declare -a win_lengths=(3000 1000)

mkdir -p logs


for win_len in "${win_lengths[@]}"; do

	echo "Win len: $win_len"
	conf_dir="conf/win_len_$win_len"
	
	# SNVs
	sbatch -J "W${win_len}_MAF0.01_snv" -o "logs/W${win_len}_MAF0.01_snv.out" -e "logs/W${win_len}_MAF0.01_snv.err" ./wgs.sh $conf_dir/config_MAF0.01_snv.yaml input_classes.txt
	sbatch -J "W${win_len}_MAF0.001_snv" -o "logs/W${win_len}_MAF0.001_snv.out" -e "logs/W${win_len}_MAF0.001_snv.err" ./wgs.sh $conf_dir/config_MAF0.001_snv.yaml input_classes.txt
	sbatch -J "W${win_len}_MAF0.0001_snv" -o "logs/W${win_len}_MAF0.0001_snv.out" -e "logs/W${win_len}_MAF0.0001_snv.err" ./wgs.sh $conf_dir/config_MAF0.0001_snv.yaml input_classes.txt


	# INDELs
	sbatch -J "W${win_len}_MAF0.01_indels" -o "logs/W${win_len}_MAF0.01_indels.out" -e "logs/W${win_len}_MAF0.01_indels.err" ./wgs.sh $conf_dir/config_MAF0.01_indels.yaml input_classes.txt
	sbatch -J "W${win_len}_MAF0.001_indels" -o "logs/W${win_len}_MAF0.001_indels.out" -e "logs/W${win_len}_MAF0.001_indels.err" ./wgs.sh $conf_dir/config_MAF0.001_indels.yaml input_classes.txt
	sbatch -J "W${win_len}_MAF0.0001_indels" -o "logs/W${win_len}_MAF0.0001_indels.out" -e "logs/W${win_len}_MAF0.0001_indels.err" ./wgs.sh $conf_dir/config_MAF0.0001_indels.yaml input_classes.txt

	# All (SNVs + INDELs)
	sbatch -J "W${win_len}_MAF0.01_all" -o "logs/W${win_len}_MAF0.01_all.out" -e "logs/W${win_len}_MAF0.01_all.err" ./wgs.sh $conf_dir/config_MAF0.01_all.yaml input_classes.txt	
	sbatch -J "W${win_len}_MAF0.001_all" -o "logs/W${win_len}_MAF0.001_all.out" -e "logs/W${win_len}_MAF0.001_all.err" ./wgs.sh $conf_dir/config_MAF0.001_all.yaml input_classes.txt	
	sbatch -J "W${win_len}_MAF0.0001_all" -o "logs/W${win_len}_MAF0.0001_all.out" -e "logs/W${win_len}_MAF0.0001_all.err" ./wgs.sh $conf_dir/config_MAF0.0001_all.yaml input_classes.txt	

done
