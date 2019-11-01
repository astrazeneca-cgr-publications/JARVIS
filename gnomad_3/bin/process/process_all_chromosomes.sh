#!/bin/bash 
#SBATCH -J process_all_chr  # Set name for current job 
#SBATCH -o out.process_all_chr  # Set job output log 
#SBATCH -e err.process_all_chr  # Set job error log 
#SBATCH --cpus-per-task=10         # Request 10 CPUs (cores) on a single node 
#SBATCH --mem=2G          # Request amount of memory 
#SBATCH -t 24:0:0            # Request 12 hours runtime


# check input arguments
if [ $# -lt 3 ]; then
	echo "[Error]: incorrect number of input arguments."
	echo ">> Expected call format: ./process_all_chr.sh [dataset: gnomad|gnomad_3|topmed] [vcf_dir: input dir with VCF files] [KEEP_PASS_ONLY: 0|1] [FILTER_SEGDUP: 0|1] [FILTER_LCR: 0|1] [population (optional): leave blank to consider all populations (otherwise: e.g. FIN, ASJ)]"
	exit
fi

dataset=$1
KEEP_PASS_ONLY=$2
FILTER_LCR=$3


version="3.0"
vcf_dir="../../vcf"
population="all"


base_out_dir="../../out"
mkdir -p $base_out_dir
out_dir=$base_out_dir/${dataset}-filtered_variant_tables


# Update output dir to store results from ./get_chr_table.sh calls
out_dir="$out_dir-$population"



if [ "$KEEP_PASS_ONLY" -eq 1 ]; then
	out_dir="$out_dir-PASS_ONLY"
fi
if [ "$FILTER_LCR" -eq 1 ]; then
	out_dir="$out_dir-NO_LCR"
fi

# Create clean output directory
rm -rf $out_dir; mkdir -p $out_dir
echo "Created out dir: $out_dir"



cnt=1
for chr in `seq 1 22`; # X Y;
do
        echo Processing chr: $chr
	input_file=$vcf_dir/gnomad.genomes.r${version}.sites.chr${chr}.vcf.bgz
 
	python -u parse_vcf.py $dataset $chr $input_file $KEEP_PASS_ONLY $FILTER_LCR $out_dir &

	if [ $cnt == 24 ]; then     
	         wait
                 cnt=1
        fi         
	cnt=$((cnt + 1))
done

wait
