#!/bin/bash

labelset_size=$1 #10000
discard_zero_values=$2


file_annot="D${labelset_size}"
if [ $discard_zero_values == "1" ]; then
	file_annot="${file_annot}.no_zeros"
fi


input_dir="most_least_conserved_files_by_class"
out_dir="${input_dir}/${file_annot}"
mkdir -p $out_dir


function label_class {

	genomic_class=$1
	echo "> $genomic_class"

	# - least conserved
	echo "Labelling least conserved regions ..."
	least_conserved_file=$out_dir/${genomic_class}.least_conserved.${file_annot}.bed
	
	# keep larger set of least-conserved regions, to sample from (as several may be in contiguous places)
	if [ $discard_zero_values == "1" ]; then
		tail_len=$(($labelset_size * 2))
		echo "Discarding regions with "0.0" conservation score ..."
			zcat < $input_dir/${genomic_class}.phastCons46way.primates.high_conf_regions.bed.sorted.gz | grep -v "0\.000000" | tail -n $tail_len | shuf -n $labelset_size | awk '{print $1"\t"$2"\t"$3"\tNon_conserved" }' > $least_conserved_file 

	else
		if [ $genomic_class == "ucne" ]; then
			tail_len=$(($labelset_size * 5))
			zcat < $input_dir/${genomic_class}.phastCons46way.primates.high_conf_regions.bed.sorted.gz | tail -n $tail_len | shuf -n $labelset_size | awk '{print $1"\t"$2"\t"$3"\tNon_conserved" }' > $least_conserved_file 
		else 
			zcat < $input_dir/${genomic_class}.phastCons46way.primates.high_conf_regions.bed.sorted.gz | grep "0\.000000" | shuf -n $labelset_size | awk '{print $1"\t"$2"\t"$3"\tNon_conserved" }' > $least_conserved_file 
		fi
	fi

	# - most conserved
	echo "Labelling most conserved regions ..."
	most_conserved_file=$out_dir/${genomic_class}.most_conserved.${file_annot}.bed

	zcat < $input_dir/${genomic_class}.phastCons46way.primates.high_conf_regions.bed.sorted.gz | grep "1\.000000" | awk '{print $1"\t"$2"\t"$3"\tConserved" }' > $most_conserved_file 
	zcat < $input_dir/${genomic_class}.phastCons46way.primates.high_conf_regions.bed.sorted.gz | grep "0\.999000" | shuf -n $labelset_size | awk '{print $1"\t"$2"\t"$3"\tConserved" }' >> $most_conserved_file

}



for genomic_class in ucne utr vista ccds intron lincrna intergenic; do
#for genomic_class in intron; do
	label_class $genomic_class &
done

wait


# Concatenate files from all classes, both most and least conserved
# - To be used for subsetting other genome-wide scores for the benchmarking
cat $out_dir/*most_conserved.${file_annot}.* > $out_dir/All_classes.conservation_annot_regions.${file_annot}.bed
cat $out_dir/*least_conserved.${file_annot}.* >> $out_dir/All_classes.conservation_annot_regions.${file_annot}.bed

cat $out_dir/All_classes.conservation_annot_regions.${file_annot}.bed | sed 's/^chr//g' > $out_dir/All_classes.conservation_annot_regions.${file_annot}.no_chr.bed
