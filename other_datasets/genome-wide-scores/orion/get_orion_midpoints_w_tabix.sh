#!/bin/bash -l
#SBATCH -J get_orion_midpoints_w_tabix  # Set name for current job
#SBATCH -o out.get_orion_midpoints_w_tabix  # Set job output log
#SBATCH -e err.get_orion_midpoints_w_tabix  # Set job error log
#SBATCH --cpus-per-task=10         # Request 9 CPUs (cores) on a single node
#SBATCH --mem=50000          # Request amount of memory
#SBATCH -t 48:0:0            # Request 24 hours runtime

orion_scores_file=plos_one_scores/all_chr/all_chr.bed.gz # orion.1001.masked.new.txt.gz
out_dir=plos-one-midpoint_bed

bed_dir="$out_dir"
bed_files=($(ls $bed_dir/*bed))



for f in "${bed_files[@]}"
do
	f=`echo $f | sed 's/.*\///'`
	echo $f
	echo "tabix $orion_scores_file -B $bed_dir/$f > $out_dir/${f}.orion.plosone.bed"
	tabix $orion_scores_file -B $bed_dir/$f > $out_dir/${f}.orion.plosone.bed &
done

wait
