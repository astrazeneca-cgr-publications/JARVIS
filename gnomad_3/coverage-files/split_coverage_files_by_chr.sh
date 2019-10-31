#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:0:0

full_cov_file="gnomad.genomes.r3.0.coverage.summary.tsv.bgz"
#full_cov_file="tmp/temp.tsv.gz"


for chr in `seq 1 22` X Y; do
	
	# ~ chr1:10001 
	echo -e "\nUnfolding converage file for chr $chr ..."

	zcat < $full_cov_file | awk -v chr="chr${chr}:" '{ if($1 ~ chr) print $0}' | sed 's/\:/\t/' | sed 's/chr//' | gzip > gnomad.genomes.r3.0.chr${chr}.coverage.txt.gz &

done

wait
