#!/bin/bash
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:0:0

#version=hg19
high_conf_regions="../../../genomic-high-confidence-regions/high_conf_genomic_regions.with_chr.bed.gz"


cnt=0
for chr in `seq 1 22`; do
	
	echo "Intersecting chr $chr with high-conf regions ..."

	intersectBed -a bed/chr${chr}.phastCons46way.primates.bed -b $high_conf_regions > bed/chr${chr}.phastCons46way.primates.high_conf_regions.bed &

	cnt=$((cnt+1))
	if [ $cnt == 1 ]; then
		cnt=0
		wait
	fi
done

wait
