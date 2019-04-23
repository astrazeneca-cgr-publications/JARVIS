#!/bin/bash
#SBATCH -J keep_gwrvis_high_conf	# Set name for current job
#SBATCH -o out.keep_gwrvis_high_conf  # Set job output log
#SBATCH -e err.keep_gwrvis_high_conf  # Set job error log
#SBATCH --cpus-per-task=4         # Request 3 CPUs (cores) on a single node
#SBATCH --mem=8G          # Request amount of memory
#SBATCH -t 24:0:0           # Request 24 hours runtime


if [ $# -ne 3 ]
then
	printf "\n>> Expected call format: [sbatch] keep_gwrvis_high_conf_regions.sh [dataset: gnomad|topmed] [input_filtered_dir] [population]\n"
	exit
fi

dataset=$1
input_filtered_dir=$2  #"gnomad-filtered_variant_tables-[all]-[filter_annot]"
population=$3

# BED file with high confidence genomic regions
include_file=../../../genomic-high-confidence-regions/high_conf_genomic_regions.bed.gz
#include_file=../../coverage-files/high_cov_bed_files-min_depth20/gnomad.whole_genome.high_cov.min_depth20.bed


filter_annot=`echo $input_filtered_dir | sed 's/.*-//' | sed 's/\///'`
#echo $filter_annot
out="../../out/${dataset}-high_conf_regions-filtered_variant_tables-${filter_annot}"
mkdir -p $out



function filter {
	i=$1
	echo Chr $i;
	
	cat $input_filtered_dir/chr${i}_${dataset}_table.all.txt.filtered | head -1 > $out/tmp0_chr${i}.txt
	# convert coordinates from VCF to 0-based for use with BED file
	cat $input_filtered_dir/chr${i}_${dataset}_table.all.txt.filtered | cut -f1 | awk '{print $1-1}' > $out/tmp1_chr${i}.txt

	paste $out/tmp1_chr${i}.txt $input_filtered_dir/chr${i}_${dataset}_table.all.txt.filtered | tail -n+2 > $out/tmp2_chr${i}.txt

	sed "s/^/${i}\t/g" $out/tmp2_chr${i}.txt > $out/tmp3_chr${i}.txt

	intersectBed -a $out/tmp3_chr${i}.txt -b $include_file > $out/tmp4_chr${i}.txt

	cat $out/tmp4_chr${i}.txt | cut --complement -f1,2 > $out/tmp5_chr${i}.txt

	cat $out/tmp0_chr${i}.txt $out/tmp5_chr${i}.txt > $out/chr${i}_${dataset}_table.${population}.txt.filtered
}

cnt=0
for k in `seq 1 22`;
do
	filter $k &

	if [ $cnt = 4 ]; then     
	         wait
                 cnt=0
        fi         
	cnt=$((cnt + 1))

	sleep 1
done
wait

# cleanup tmp files
rm -rf $out/tmp*
