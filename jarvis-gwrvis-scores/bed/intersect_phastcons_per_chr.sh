#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=60G 
#SBATCH --time=24:0:0

function subset_chr {

	chr=$1
	echo "Chr: $chr started..."
	intersectBed -wo -a JARVIS.prediction_scores.no_genomic_classes.bed -b ../../other_datasets/conservation/phastCons46way_primates/bed/chr${chr}.phastCons46way.primates.high_conf_regions.bed.gz | gzip > jarvis-phastcons-output/JARVIS_phastcons.full_table.chr${chr}.bed.gz

	echo "Chr: $chr [Done]"

}


mkdir -p jarvis-phastcons-output

cnt=0
for chr in `seq 1 22`; do

	subset_chr $chr &

	cnt=$((cnt+1))
	if [ "$cnt" == 2 ]; then
		cnt=0
		wait
	fi
	
done

wait
