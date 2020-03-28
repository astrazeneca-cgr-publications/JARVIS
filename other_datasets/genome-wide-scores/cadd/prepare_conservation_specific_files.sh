#!/bin/bash
#SBATCH -t 24:0:0
#SBATCH -n 1
#SBATCH --mem=4G

labelset_size=$1
discard_zero_values=$2

file_annot="D${labelset_size}"
if [ $discard_zero_values == "1" ]; then
	file_annot="${file_annot}.no_zeros"
fi

conservation_ref_file="../../conservation/phastCons46way_primates/most_least_conserved_files_by_class/${file_annot}/All_classes.conservation_annot_regions.${file_annot}.no_chr.bed"

full_score_ref_file="cadd.whole_genome.all_raw_scores.vcf.bgz"
score="cadd"



# > Pathogenic variant-specific files
tabix $full_score_ref_file -B $conservation_ref_file | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > ${score}.conservation_annotation.${file_annot}.bed


