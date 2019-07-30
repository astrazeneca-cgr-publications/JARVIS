#!/bin/sh

out_dir=$1
score=$2


out_dir="../"$out_dir"/full_genome_out"
benchmark_dir=$out_dir"/scores_benchmarking"
mkdir -p $benchmark_dir
benchmark_dir="$benchmark_dir/$score"
mkdir -p $benchmark_dir


intersectBed -a $out_dir/BED/full_genome.All_genomic_classes.bed -b ../../other_datasets/clinvar/clinvar.pathogenic.bed > $benchmark_dir/$score.clinvar_pathogenic.bed
intersectBed -a $out_dir/BED/full_genome.All_genomic_classes.bed -b ../../other_datasets/clinvar/clinvar.benign.bed > $benchmark_dir/$score.clinvar_benign.bed

echo "python run_score_benchmarking.py $benchmark_dir/$score.clinvar_pathogenic.bed $benchmark_dir/$score.clinvar_benign.bed"
python run_score_benchmarking.py $benchmark_dir/$score.clinvar_pathogenic.bed $benchmark_dir/$score.clinvar_benign.bed

exit


non_cod_include_classes="intergenic intron lincrna mature_mirna ucne utr vista"
coding_include_classes="omim-HI tolerant intolerant ccds"

for cl in $non_cod_include_classes $coding_include_classes; do
	echo $cl

	pathogenic_scores_file="${cl}.clinvar_pathogenic.${score}.txt"
	benign_scores_file="${cl}.clinvar_benign.${score}.txt"

	# Subset ClinVar pathogenic/benign variants per genomic class
	cat $benchmark_dir/$score.clinvar_pathogenic.bed | grep -E $cl > $pathogenic_scores_file
	cat $benchmark_dir/$score.clinvar_benign.bed | grep -E $cl > $benign_scores_file

	echo $pathogenic_scores_file
	exit

	# get BED files per genomic class
	cat $benchmark_dir/$score.clinvar_pathogenic.bed | grep -E $cl | cut --complement -f4 | sed 's/chr//' > ${cl}.clinvar_pathogenic.bed
	cat  | grep -E $cl | cut --complement -f4 | sed 's/chr//' > ${cl}.clinvar_benign.bed

	#python plot_clinvar_densities.py $pathogenic_scores_file $benign_scores_file $cl

	echo "python run_logistic_regression.py $pathogenic_scores_file $benign_scores_file $cl"
	python run_logistic_regression.py $pathogenic_scores_file $benign_scores_file $cl


done
