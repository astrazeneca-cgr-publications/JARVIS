#!/bin/bash -l
#SBATCH -J keep_gwrvis_high_conf	# Set name for current job
#SBATCH -o out.keep_gwrvis_high_conf  # Set job output log
#SBATCH -e err.keep_gwrvis_high_conf  # Set job error log
#SBATCH --cpus-per-task=8         # Request 1 CPU (core) on a single node
#SBATCH --mem=30000          # Request amount of memory
#SBATCH -t 48:0:0           # Request 48 hours runtime

input_dir=$1
dataset=$2	# gnomad or topmed

if [ $# -ne 2 ]
then
	echo "Please provide:"
	echo "    - an 'input_dir' arg, containing VCF files with (filtered) variants per chromosome and"
	echo "    - a 'dataset' arg (gnomad or topmed) and"
	printf "\n>> Expected call format: [sbatch] keep_gwrvis_high_confidence_regions.sh dataset(gnomad|topmed)\n"
	exit
fi


if [[ "$dataset" == gnomad || "$dataset" == topmed ]]
then
echo	"[OK] Dataset: ${dataset}"
else
	echo "Please provide a valid value for the 'dataset' argument."
	echo "    - Accepted values: 'gnomad' or 'topmed'"
	printf "\n>> Expected call format: [sbatch] keep_gwrvis_high_confidence_regions.sh dataset(gnomad|topmed)\n"
	exit
fi


mkdir out
out="out/${dataset}-high_confidence-${input_dir}"
include_file=../genomic-high-confidence-regions/high_conf_genomic_regions.bed.gz 

mkdir -p $out

function filter {
	i=$1

	echo Chr $i;
	
	cat $input_dir/chr${i}_${dataset}_table.all.txt.filtered | head -1 > $out/tmp_chr${i}.header


	# convert coordinates from VCF to 0-based for use with BED file
	tail -n+2 $input_dir/chr${i}_${dataset}_table.all.txt.filtered | sed 's/^chr//' | awk -v chrom="$i" '{print chrom"\t"$1-1"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8 }' > $out/tmp_chr${i}.bed
	

	# intersect with high-confidence genomic regions
	intersectBed -a $out/tmp_chr${i}.bed -b $include_file | sed 's/^/chr/' > $out/tmp_chr${i}.high_conf.bed

	# convert BED back to VCF
	cat $out/tmp_chr${i}.high_conf.bed | cut --complement -f1,2 > $out/tmp_chr${i}.high_conf.vcf

	# Save original header and filtered VCF to a file
	cat $out/tmp_chr${i}.header $out/tmp_chr${i}.high_conf.vcf > $out/chr${i}_${dataset}_table.all.txt.filtered
}

cnt=0
for k in `seq 1 22`;
do
	filter $k &

	if [ $cnt = 7 ]; then     
	         wait
                 cnt=0
        fi         
	cnt=$((cnt + 1))

	sleep 1
done

wait

rm -rf $out/tmp*
