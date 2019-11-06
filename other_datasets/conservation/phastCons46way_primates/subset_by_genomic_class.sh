#!/bin/bash
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=1
#SBATCH --time=24:0:0

out_dir="bed-by_class"
mkdir -p $out_dir


function subset_for_class {

	genomic_class=$1

	genomic_class_bed_dir="genomic_classes_mutually_excl_bed-topmed"
	cur_genomic_class_bed="${genomic_class_bed_dir}/gwrvis.${genomic_class}.mutually_excl.bed"

	out_file="$out_dir/${genomic_class}.phastCons46way.primates.high_conf_regions.bed"
	rm -f $out_file

	for chr in `seq 1 22`; do
		echo "Chr: $chr"

		cur_chr_phastcons_bed="bed/chr${chr}.phastCons46way.primates.high_conf_regions.bed.gz"

		tabix $cur_chr_phastcons_bed -B $cur_genomic_class_bed >> $out_file
	done

	cat $out_file | sort -k4,4rn | gzip > ${out_file}.sorted.gz
	#rm $out_file
}



cnt=0
for genomic_class in ucne utr vista ccds intron lincrna intergenic; do

	echo "> $genomic_class"
	subset_for_class $genomic_class &

	cnt=$((cnt+1))
	if [ $cnt == 3 ]; then
		cnt=0
		wait
	fi
done

wait
