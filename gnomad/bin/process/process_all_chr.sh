#!/bin/bash -l
#SBATCH -J process_all_chr  # Set name for current job 
#SBATCH -o out.process_all_chr  # Set job output log 
#SBATCH -e err.process_all_chr  # Set job error log 
#SBATCH --cpus-per-task=10         # Request 10 CPUs (cores) on a single node 
#SBATCH --mem=4G          # Request amount of memory 
#SBATCH -t 12:0:0            # Request 12 hours runtime


# check input arguments
if [ $# -lt 4 ]; then
	echo "[Error]: incorrect number of input arguments."
	echo ">> Expected call format: ./process_all_chr.sh [input_dir_with_vcf_files] [KEEP_PASS_ONLY: 0|1] [FILTER_SEGDUP: 0|1] [FILTER_LCR: 0|1] [population (optional): leave blank to consider all populations (otherwise: e.g. FIN, ASJ)]"
	exit
fi

vcf_dir=$1
KEEP_PASS_ONLY=$2 
FILTER_SEGDUP=$3
FILTER_LCR=$4
population=$5


filter_by_subpopulation()
{
	population=$1
	
	in_vcf_dir=$vcf_dir
	out_vcf_dir="$vcf_dir-$population"
	rm -rf $out_vcf_dir/*; mkdir -p $out_vcf_dir;

	# Also update $vcf_dir to be passed onto ./get_chr_table.sh
	vcf_dir="$vcf_dir-$population"

	cnt=1
	# BETA
	for i in `seq 1 22`;
	do
		echo Processing chr: $i
		python subpop_vcf_filter.py $i $population $in_vcf_dir $out_vcf_dir &

		if [ $cnt = 10 ]; then     
			 wait
			 cnt=1
		fi         
		cnt=$((cnt + 1))
	done

	i=X
	echo Processing chr: $i
	python subpop_vcf_filter.py $i $population $in_vcf_dir $out_vcf_dir &

	wait

}

# Filter by subpopulation (if applicable) and create respective output dir
out_dir=filtered_variant_tables
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


#for file in $vcf_dir/*.vcf; do
#	echo "Zipping file: $file ...";
#	gzip -S ".bgz" $file &
#done
#wait


cnt=1
for i in `seq 1 22`;
do
        echo Processing chr: $i
	./get_chr_table.sh $vcf_dir $i $KEEP_PASS_ONLY $FILTER_SEGDUP $FILTER_LCR $out_dir $population &

	if [ $cnt = 9 ]; then     
	         wait
                 cnt=1
        fi         
	cnt=$((cnt + 1))
done

i=X
echo Processing chr: $i
./get_chr_table.sh $vcf_dir $i $KEEP_PASS_ONLY $FILTER_SEGDUP $FILTER_LCR $out_dir $population &
wait
