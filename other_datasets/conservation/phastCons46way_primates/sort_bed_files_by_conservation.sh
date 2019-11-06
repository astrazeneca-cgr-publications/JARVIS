#!/bin/sh
#SBATCH -t 24:0:0
#SBATCH -n 2
#SBATCH --mem-per-cpu=8G

cnt=0
for chr in `seq 1 22`; do
	
	echo "Chr: $chr"

	cat bed/chr${chr}.phastCons46way.primates.bed | sort -k4,4n > bed/chr${chr}.phastCons46way.primates.sorted.bed &

	# Set max num of threads
	if [ $cnt == 4 ]; then
		cnt=0
		wait
	fi
	cnt=$((cnt+1))

done

wait
