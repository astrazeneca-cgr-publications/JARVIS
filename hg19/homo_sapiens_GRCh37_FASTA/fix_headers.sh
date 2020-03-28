#!/bin/bash

version="37.75"

for chr in `seq 1 22` X Y; do
	echo "Fixing header for chr $chr ..."
 
	zcat < Homo_sapiens.GRCh${version}.dna.chromosome.${chr}.fa.gz | awk -v chr="$chr" '{ if($0 ~ /^>/ ) print ">chr"chr; else print $0 }' | gzip > Homo_sapiens.GRCh${version}.dna.chromosome.${chr}.fa.gz.tmp & 

done

wait

	
echo "Renaming files with fixed headers back to the original ones ..."
for chr in `seq 1 22` X Y; do

	mv Homo_sapiens.GRCh${version}.dna.chromosome.${chr}.fa.gz.tmp Homo_sapiens.GRCh${version}.dna.chromosome.${chr}.fa.gz 

done

echo "Done."
