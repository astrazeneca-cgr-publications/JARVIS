#!/bin/bash
#SBATCH -t 24:0:0
#SBATCH -n 1
#SBATCH --mem=4G


#declare -a pathogenic_resources=(clinvar hgmd)
#declare -a pathogenic_resources=(ncER-training ncER-gwas_catalog ncER-test ncER-mendelian ncER-generalization_other ncER-generalization_ncRNA)
#declare -a pathogenic_resources=(ncER-training-full ncER-training-NoChr5_6_8 ncER-validation-Chr5_6_8)
#declare -a benign_resources=(clinvar denovodb denovodb_nonSSC topmed_uniq)
#declare -a pathogenic_resources=(clinvar-test-Chr568 clinvar-training-NoChr568)
declare -a pathogenic_resources=(ncER-mendelian-clean ncER-generalization_ncRNA-clean ncER-generalization_other-clean ncER-test-clean ncER-gwas_catalog-clean)

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
