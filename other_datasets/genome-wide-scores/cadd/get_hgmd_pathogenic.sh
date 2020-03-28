#!/bin/bash
#SBATCH -t 24:00:00

tabix cadd.whole_genome.all_raw_scores.vcf.bgz -B ../../hgmd/hgmd.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > cadd.hgmd_pathogenic.bed   

cat cadd.hgmd_pathogenic.bed cadd.clinvar_benign.bed > cadd.hgmd.bed
