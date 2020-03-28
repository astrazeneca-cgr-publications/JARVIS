#!/bin/sh

phenotypes=(`cat denovo-db.non-ssc-samples.variants.bed | cut -f5 | sort | uniq`)
 

for pheno in "${phenotypes[@]}"; do 
	echo $pheno":"

	cat denovo-db.non-ssc-samples.variants.bed | grep $pheno | wc -l
done
