#!/bin/bash

module load module load Java/1.8.0_181

version="37.75"
echo "gatk CreateSequenceDictionary ..."
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh${version}.dna.all_chromosome.fa


# create .fai index file
echo "samtools faidx ..."
samtools faidx Homo_sapiens.GRCh${version}.dna.all_chromosome.fa


# Retrieve sequences from indexed human genome
v=`echo $version | sed 's/\..*//'`
echo "faTwoBit ..."
faToTwoBit Homo_sapiens.GRCh${v}.dna.all_chromosome.fa hsa${v}.2bit

echo "twoBitInfo ..."
twoBitInfo hsa${v}.2bit stdout | sort -k2rn > hsa${v}.chrom.sizes

# Test: 
twoBitToFa hsa${v}.2bit:1:1-10 /dev/stdout
