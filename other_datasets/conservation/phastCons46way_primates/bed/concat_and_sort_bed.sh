#!/bin/bash
#SBATCH --mem=40G
#SBATCH -t 24:00:0


final_file="All_chr.phastCons46way.primates.sorted.bed"

#for chr in `seq 1 22`; do
#	cat chr${chr}.phastCons46way.primates.bed >> $final_file
#done

echo "Zipping All_chr file..."
zcat < tmp.${final_file}.gz | bgzip > ${final_file}.gz

"Building tabix index..."
tabix -p bed ${final_file}.gz
