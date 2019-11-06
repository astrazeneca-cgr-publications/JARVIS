#!/bin/bash -l 
#SBATCH -J gatk_liftover  # Set name for current job 
#SBATCH -o out.gatk_liftover  # Set job output log 
#SBATCH -e err.gatk_liftover  # Set job error log 
#SBATCH --cpus-per-task=2         # Request 1 CPU (core) on a single node 
#SBATCH --mem=40000          # Request amount of memory 
#SBATCH -t 48:0:0           # Request 48 hours runtime

vcf_dir=vcf

gatk LiftoverVcf -C $vcf_dir/hg38ToHg19.over.chain -I $vcf_dir/bravo-dbsnp-all.vcf -O $vcf_dir/bravo-dbsnp-all.hg19.vcf -R ../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.fa --REJECT=${vcf_dir}/rejected_topmed_entries_from_liftover.vcf --WARN_ON_MISSING_CONTIG=true --LIFTOVER_MIN_MATCH=0.90
