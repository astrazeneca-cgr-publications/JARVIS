#!/bin/bash
#SBATCH -J sort_big_bed.faster_2
#SBATCH -o sort_big_bed.faster_2.out
#SBATCH -e sort_big_bed.faster_2.err
#SBATCH --cpus-per-task=35
#SBATCH --mem=200G
#SBATCH -t 24:00:00


config_file=$1

out_dir=`python custom_utils.py $config_file`
unsorted_file="${out_dir}/full_gwrvis_single_nt.unsorted.bed"
echo "Unsorted file: $unsorted_file"

sorted_file="${unsorted_file}.sorted.faster_2.bed"


sort --parallel=100 -k1,1V -k2,2n -o $sorted_file $unsorted_file 


echo "Sorted file: $sorted_file"
