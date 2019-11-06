#!/bin/sh
#SBATCH --time=48:0:0

#tabix -p vcf DANN_whole_genome_SNVs.tsv.bgz

tabix DANN_whole_genome_SNVs.tsv.bgz -B ../../clinvar/clinvar.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > dann.clinvar_pathogenic.bed & 
tabix DANN_whole_genome_SNVs.tsv.bgz -B ../../clinvar/clinvar.benign.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > dann.clinvar_benign.bed

wait
cat dann.clinvar_pathogenic.bed dann.clinvar_benign.bed > dann.clinvar.bed
