#!/bin/sh
#SBATCH --time=48:0:0

zcat < whole_genome_SNVs.tsv.gz | grep -v '#' | cut -f1,2,3,4,5 | bgzip -c > cadd.whole_genome.all_raw_scores.vcf.bgz

tabix -p vcf cadd.whole_genome.all_raw_scores.vcf.bgz
