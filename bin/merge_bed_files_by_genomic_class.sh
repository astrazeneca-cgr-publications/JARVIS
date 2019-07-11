#!/bin/sh

out_dir=$1
input_classes_file=$2
readarray -t input_classes < $input_classes_file

gwrvis_distr_dir="${out_dir}/gwrvis_distribution"


for class in "${input_classes[@]}"; do
	echo $cl
	genomic_class_out_file="${gwrvis_distr_dir}/gwrvis.${class}.mutually_excl.bed"

	echo "cat ${gwrvis_distr_dir}/BED/gwrvis_scores_chr*.genomic_coords.${class}.bed > $genomic_class_out_file"
	cat ${gwrvis_distr_dir}/BED/gwrvis_scores_chr*.genomic_coords.${class}.bed > $genomic_class_out_file

	class_total_nt_lenth=`cat $genomic_class_out_file | awk '{sum+=$3-$2} END {print sum}'`
	echo -e "$class total nt length: $class_total_nt_lenth \n"
done



