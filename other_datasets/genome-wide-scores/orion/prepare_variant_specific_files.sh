#!/bin/bash
#SBATCH -t 24:0:0
#SBATCH -n 1
#SBATCH --mem=4G


#declare -a pathogenic_resources=(clinvar hgmd)
#declare -a pathogenic_resources=(ncER-training ncER-gwas_catalog ncER-test ncER-mendelian ncER-generalization_other ncER-generalization_ncRNA)
declare -a pathogenic_resources=(ncER-training-full)
#declare -a benign_resources=(clinvar denovodb denovodb_nonSSC topmed_uniq)

full_score_ref_file="orion.1001.masked.new.txt.gz"
score="orion"



# > Pathogenic variant-specific files
for resource in "${pathogenic_resources[@]}"; do
	echo "pathogenic: $resource"
	echo ../../variant_annotation-hg19/${resource}/${resource}.pathogenic.no_chr_prefix.bed

	tabix $full_score_ref_file -B ../../variant_annotation-hg19/${resource}/${resource}.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > ${score}.${resource}_pathogenic.bed &
done


# > Benign variant-specific files
for resource in "${benign_resources[@]}"; do
	echo "benign: $resource"
	tabix $full_score_ref_file -B ../../variant_annotation-hg19/${resource}/${resource}.benign.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > ${score}.${resource}_benign.bed &
done

wait
