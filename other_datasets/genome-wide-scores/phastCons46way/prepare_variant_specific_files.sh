#!/bin/bash
#SBATCH -t 48:0:0
#SBATCH -n 4
#SBATCH --mem=4G

declare -a pathogenic_resources=(clinvar hgmd)
declare -a benign_resources=(clinvar denovodb denovodb_nonSSC topmed_uniq)

full_score_ref_file="All_chromosomes.phastCons46way.bed.bgz"
score="phastCons46way"


# > Pathogenic variant-specific files
for resource in "${pathogenic_resources[@]}"; do
	tabix $full_score_ref_file -B ../../variant_annotation/${resource}/${resource}.pathogenic.bed > ${score}.${resource}_pathogenic.bed &
done


# > Benign variant-specific files
for resource in "${benign_resources[@]}"; do
	tabix $full_score_ref_file -B ../../variant_annotation/${resource}/${resource}.benign.bed > ${score}.${resource}_benign.bed &
done


wait
