#!/bin/bash
#SBATCH -t 48:0:0
#SBATCH -n 1
#SBATCH --mem=4G

declare -a pathogenic_resources=(ncER-training ncER-gwas_catalog ncER-test ncER-mendelian ncER-generalization_other ncER-generalization_ncRNA)
#declare -a pathogenic_resources=(clinvar hgmd)
#declare -a benign_resources=(clinvar denovodb denovodb_nonSSC topmed_uniq)

full_score_ref_file="All_chromosomes.ncRVIS.bed"
score="ncRVIS"


# > Pathogenic variant-specific files
for resource in "${pathogenic_resources[@]}"; do
	intersectBed -a $full_score_ref_file -b ../../variant_annotation-hg19/${resource}/${resource}.pathogenic.bed > ${score}.${resource}_pathogenic.bed &
done


# > Benign variant-specific files
for resource in "${benign_resources[@]}"; do
	intersectBed -a $full_score_ref_file -b ../../variant_annotation-hg19/${resource}/${resource}.benign.bed > ${score}.${resource}_benign.bed &
done


wait
