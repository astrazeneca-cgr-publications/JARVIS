#!/bin/bash 
#SBATCH -J parse_all_chromosomes	# Set name for current job 
#SBATCH -o out.parse_all_chromosomes	# Set job output log 
#SBATCH -e err.parse_all_chromosomes	# Set job error log 
#SBATCH --cpus-per-task=1         # Request 22 CPUs (cores) on a single node 
#SBATCH --mem=4000          # Request amount of memory 
#SBATCH -t 24:0:0            # Request 24 hours runtime

config_log=$1

cnt=0
#for i in `seq 1 22`;
for i in `seq 21 21`;
do
	echo Running parse-job for chr: $i
	python filtered_chr_parser.py $i $config_log #&

	if [ $cnt = 8 ]; then
		wait
		cnt=0
	fi
	cnt=$((cnt + 1))
done

#i=X
#echo Running parse-job for chr: $i
#python full_genome_filtered_chr_parser.py $i $config_log
wait
