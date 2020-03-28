#!/bin/bash 

version="37.75"

rm -f Homo_sapiens.GRCh${version}.dna.all_chromosome.fa

for i in `seq 1 22` X Y;
do

	cat Homo_sapiens.GRCh${version}.dna.chromosome.${i}.fa.gz >> Homo_sapiens.GRCh${version}.dna.all_chromosome.fa

done


mv Homo_sapiens.GRCh${version}.dna.all_chromosome.fa Homo_sapiens.GRCh${version}.dna.all_chromosome.fa.gz


# unzip concatenated fa.gz
gunzip Homo_sapiens.GRCh${version}.dna.all_chromosome.fa.gz

cat Homo_sapiens.GRCh${version}.dna.all_chromosome.fa | sed 's/chr//g' > Homo_sapiens.GRCh${version}.dna.all_chromosome.fa.tmp
mv Homo_sapiens.GRCh${version}.dna.all_chromosome.fa.tmp Homo_sapiens.GRCh${version}.dna.all_chromosome.fa
