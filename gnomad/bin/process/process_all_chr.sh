#!/bin/bash 
#SBATCH -J process_all_chr  # Set name for current job 
#SBATCH -o out.process_all_chr  # Set job output log 
#SBATCH -e err.process_all_chr  # Set job error log 
#SBATCH --cpus-per-task=25         # Request 10 CPUs (cores) on a single node 
#SBATCH --mem=4G          # Request amount of memory 
#SBATCH -t 24:0:0            # Request 12 hours runtime


# check input arguments
if [ $# -lt 5 ]; then
	echo "[Error]: incorrect number of input arguments."
	echo ">> Expected call format: ./process_all_chr.sh [dataset: gnomad|topmed] [vcf_dir: input dir with VCF files] [KEEP_PASS_ONLY: 0|1] [FILTER_SEGDUP: 0|1] [FILTER_LCR: 0|1] [population (optional): leave blank to consider all populations (otherwise: e.g. FIN, ASJ)]"
	exit
fi

dataset=$1
vcf_dir=$2
KEEP_PASS_ONLY=$3 
FILTER_SEGDUP=$4
FILTER_LCR=$5
population=$6

base_out_dir="../../out"
mkdir -p $base_out_dir
out_dir=$base_out_dir/${dataset}-filtered_variant_tables


filter_by_subpopulation()
{
	population=$1
	
	in_vcf_dir=$vcf_dir
	out_vcf_dir="$vcf_dir-$population"
	rm -rf $out_vcf_dir/*; mkdir -p $out_vcf_dir;

	# Also update $vcf_dir to be passed onto ./get_chr_table.sh
	vcf_dir="$vcf_dir-$population"

	cnt=1
	#for i in `seq 24 24`; # For debugging (24: artefactual chr of small size)
	for i in `seq 1 22`;
	do
		echo Processing chr: $i
		echo python subpop_vcf_filter.py $i $population $in_vcf_dir $out_vcf_dir
		python subpop_vcf_filter.py $i $population $in_vcf_dir $out_vcf_dir &

		if [ $cnt = 10 ]; then     
			 wait
			 cnt=1
		fi         
		cnt=$((cnt + 1))
	done

	#i=X
	#echo Processing chr: $i
	#python subpop_vcf_filter.py $i $population $in_vcf_dir $out_vcf_dir &
	wait

}

# Filter by subpopulation (if applicable) and create respective output dir
if [ -n "$population" ]; then

	echo "Filtering VCF by population: $population"
	# Filter VCF files based on selected sub-population 
	filter_by_subpopulation $population
else
	population='all'
fi


# Update output dir to store results from ./get_chr_table.sh calls
out_dir="$out_dir-$population"


if [ "$KEEP_PASS_ONLY" -eq 1 ]; then
	out_dir="$out_dir-PASS_ONLY"
fi
if [ "$FILTER_SEGDUP" -eq 1 ]; then
	out_dir="$out_dir-NO_SEGDUP"
fi
if [ "$FILTER_LCR" -eq 1 ]; then
	out_dir="$out_dir-NO_LCR"
fi

# Create clean output directory
rm -rf $out_dir; mkdir -p $out_dir
echo "Created out dir: $out_dir"




cnt=1
for i in `seq 1 22`;
do
        echo Processing chr: $i
	./get_chr_table.sh $vcf_dir $i $KEEP_PASS_ONLY $FILTER_SEGDUP $FILTER_LCR $out_dir $population &

	if [ $cnt = 25 ]; then     
	         wait
                 cnt=1
        fi         
	cnt=$((cnt + 1))
done

#i=X # Not available in r2.1.1
#echo Processing chr: $i
#./get_chr_table.sh $vcf_dir $i $KEEP_PASS_ONLY $FILTER_SEGDUP $FILTER_LCR $out_dir $population &
wait
