#!/bin/bash
#SBATCH --mem=100G    #250G for CADD
#SBATCH -n 4
#SBATCH -t 24:0:0

score=$1


function main () {

	score=$1
	sv_class=$2
	python process_scores_per_sv.py $score $sv_class

}



sv_classes=('utr' 'intergenic' 'lof' 'promoter' 'copy_gain' 'dup_partial' 'inv_span' 'dup_lof' 'intronic')

for cl in "${sv_classes[@]}"; do

	echo "Class: $cl"
	main $score $cl 

done
