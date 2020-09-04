#!/bin/bash
#SBATCH -J build_tabix 
#SBATCH -o build_tabix.out
#SBATCH -e build_tabix.err
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH -t 24:00:00


config_file=$1

out_dir=`python custom_utils.py $config_file`
echo $out_dir



single_nt_gwrvis_per_chr_dir="$out_dir/single_nt_gwrvis_per_chr"




function build_tabix_per_chr {
	
	chr=$1

	cur_file=${single_nt_gwrvis_per_chr_dir}/gwrvis_single_nt.chr${chr}.bed
	#echo $cur_file

	bgzip -c $cur_file > ${cur_file}.gz

	tabix -p bed ${cur_file}.gz
	
	echo "Done - chr$chr"
}



for chr in `seq 1 22`; do
	echo "chr$chr"

	build_tabix_per_chr $chr $in_file &
done
wait





