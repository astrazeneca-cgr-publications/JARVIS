#!/bin/bash 
#SBATCH -J parse_all_chromosomes	# Set name for current job 
#SBATCH -o out.parse_all_chromosomes	# Set job output log 
#SBATCH -e err.parse_all_chromosomes	# Set job error log 
#SBATCH --cpus-per-task=11         # Request 8 CPUs (cores) on a single node 
#SBATCH --mem=14G          # Request amount of memory 
#SBATCH -t 24:0:0            # Request 24 hours runtime

module load libpng/1.6.23-foss-2017a

config_log=$1

cnt=0
for i in `seq 1 22`;
do
	echo Running parse-job for chr: $i
	python -u filtered_chr_feature_extractor.py $i $config_log &

	if [ $cnt = 10 ]; then
		wait
		cnt=0
	fi
	cnt=$((cnt + 1))

	sleep 1
done

#i=X
#echo Running parse-job for chr: $i
#python full_genome_filtered_chr_parser.py $i $config_log
wait
