#!/bin/bash
#SBATCH -J compile_singlent_gwrvis
#SBATCH -o W5000.compile_singlent_gwrvis.out
#SBATCH -e W5000.compile_singlent_gwrvis.err
#SBATCH --cpus-per-task=22
#SBATCH --mem=100G
#SBATCH -t 24:00:00


config_file=$1

out_dir=`python custom_utils.py $config_file`
echo $out_dir


# calc offset limits
win_len=`cat $config_file | grep "win_len:" | sed 's/^[ \t]*win_len: //' | sed 's/[ ]*$//'`
left_offset=-$((win_len/2))
right_offset=$((win_len/2)) 
right_offset=$((right_offset-1))


single_nt_gwrvis_per_chr_dir="$out_dir/single_nt_gwrvis_per_chr"
rm -rf $single_nt_gwrvis_per_chr_dir
mkdir $single_nt_gwrvis_per_chr_dir







function get_gwrvis_per_offset_per_chr {
	
	chr=$1
	in_file=$2

	out_file=${single_nt_gwrvis_per_chr_dir}/gwrvis_single_nt.chr${chr}.bed
	#echo $out_file
	

	# subset for particular chr
	tail -n+2 $in_file | cut --complement -f1 | awk -v chrom="$chr" '{if($1 == "chr"chrom) print $0}' >> ${out_file}.cur_output


	# append to overall chromosome-file	
	if [ ! -f $out_file ]; then
		mv ${out_file}.cur_output  $out_file
	else
		paste -d = $out_file ${out_file}.cur_output > ${out_file}.tmp
	
		mv ${out_file}.tmp $out_file
		rm ${out_file}.cur_output
	fi
}



cnt=0
for single_nt_offset in $( seq $left_offset $right_offset); do

	echo "single-nt offset: $single_nt_offset"

	cur_dir="$out_dir/gwrvis_scores/single_nt_offset_${single_nt_offset}"
	in_file=${cur_dir}/full_genome.all_gwrvis.single_nt_offset_${single_nt_offset}.bed
	#echo $in_file


	for chr in `seq 1 22`; do
		get_gwrvis_per_offset_per_chr $chr $in_file &
	done
	wait


	#if [ $cnt == 10 ]; then
	#	break
	#fi
	#cnt=$((cnt+1))

done




function expand_chr_file_by_delim {

	cur_out_file=${single_nt_gwrvis_per_chr_dir}/gwrvis_single_nt.chr${chr}.bed

	sed 's/=/\n/g' $cur_out_file > ${cur_out_file}.tmp
	
	mv ${cur_out_file}.tmp $cur_out_file

	echo -e "Final chr$chr file:\n$cur_out_file"

}

for chr in `seq 1 22`; do
	echo "Expanding sorted out file for chr: $chr"
	expand_chr_file_by_delim $chr &
done
wait
