#!/bin/sh
#SBATCH -t 24:00:00
#SBATCH -n 4
#SBATCH --mem=4G

declare -a scores=('gwrvis') #('jarvis' 'ncER_10bp' 'linsight' 'orion')
declare -a pathogenic_classes=('intergenic' 'lof' 'dup_lof' 'inv_span' 'dup_partial' 'utr' 'copy_gain' 'promoter' 'intronic')


function convert_bed {

	score=$1
	class=$2
	printf "\nScore: $score - Class: $class\n"

	base_dir="${score}-sv-bed"
	cur_file="$base_dir/${score}.SV.${class}.with_coords.bed"
	out_file=`echo $cur_file | sed 's/bed$/1_based/'`

	cat $cur_file | cut -f1,3,4 | sed 's/chr//g' > $out_file

	
	printf "\nScore: $score - Class: $class -- DONE\n"
}


for score in "${scores[@]}" ; do

	for class in "${pathogenic_classes[@]}"; do
		convert_bed $score $class &
	done		
done

wait		
