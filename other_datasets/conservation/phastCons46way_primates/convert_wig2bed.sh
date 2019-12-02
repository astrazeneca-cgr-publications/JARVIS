#!/bin/bash
#SBATCH -t 24:0:0
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

mkdir -p bed 

cnt=0
for chr in `seq 1 22`; do
	echo "Chr $chr"

	zcat < wigFix/chr${chr}.phastCons46way.primates.wigFix.gz | wig2bed - | cut -f1,2,3,5 > bed/chr${chr}.phastCons46way.primates.bed &

	if [ $cnt == 7 ]; then
		cnt=0
		wait
	fi

	cnt=$((cnt+1))
done

wait
