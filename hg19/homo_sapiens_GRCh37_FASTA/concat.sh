#!/bin/bash 

rm -f Homo_sapiens.GRCh37.75.dna.all_chromosome.fa

for i in `seq 1 22`;
do
	cat Homo_sapiens.GRCh37.75.dna.chromosome.${i}.fa.gz >> Homo_sapiens.GRCh37.75.dna.all_chromosome.fa
done

i=X
cat Homo_sapiens.GRCh37.75.dna.chromosome.${i}.fa.gz >> Homo_sapiens.GRCh37.75.dna.all_chromosome.fa

i=Y
cat Homo_sapiens.GRCh37.75.dna.chromosome.${i}.fa.gz >> Homo_sapiens.GRCh37.75.dna.all_chromosome.fa

mv Homo_sapiens.GRCh37.75.dna.all_chromosome.fa Homo_sapiens.GRCh37.75.dna.all_chromosome.fa.gz
