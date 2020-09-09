#!/bin/sh
#SBATCH -o out.split_1_based
#SBATCH -t 10:00:00
#SBATCH -n 20
#SBATCH --mem=4G

score=$1

echo "Score: $score"
score_base_dir="${score}-sv-bed"


sv_classes=('intergenic' 'utr' 'lof' 'promoter' 'copy_gain' 'dup_partial' 'inv_span' 'dup_lof' 'intronic')
#sv_classes=('inv_span' 'dup_lof' 'intronic')


for class in "${sv_classes[@]}"; do
	echo "Class: $class"
	out_dir="${score_base_dir}/${class}-1_based"
	mkdir -p $out_dir

	in_file="${score_base_dir}/${score}.SV.${class}.with_coords.1_based"

	for chr in `seq 1 22`; do
		out_file="${out_dir}/${score}.SV.${class}.with_coords.chr${chr}.1_based"
		cat $in_file | awk -v chr="$chr" '{if($1 == chr) print $2"\t"$3}' > $out_file &
	done
done
wait

