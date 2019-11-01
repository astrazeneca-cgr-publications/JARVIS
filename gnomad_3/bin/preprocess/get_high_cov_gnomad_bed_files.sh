#!/bin/bash
#SBATCH -J 20mindepth_get_high_cov_bed
#SBATCH -o out.get_high_cov_bed_mindepth_20
#SBATCH -e err.get_high_cov_bed_mindepth_20
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00


# ------------ Read input arguments ------------
if [ $# -eq 0 ]; then
	echo "[Error]: No arguments supplied" 
	echo ">>> Please provide a min_depth (coverage) value"
	exit
fi

min_depth=$1	# 20%-percentile for gnomad mean coverage is 29.87, so suggested default min_depth: 30 -- However, x20 is already sufficient for WGS data so using this as my current threshold.
# ----------------------------------------------


version="3.0"

# ------------ Define output dirs & files ---------------
cov_file_dir="../../coverage-files"
out_dir="${cov_file_dir}/high_cov_bed_files-min_depth$min_depth"
mkdir -p $out_dir

out_all_chrom_bed_file="$out_dir/gnomad.whole_genome.high_cov.min_depth${min_depth}.bed"
rm -f $out_all_chrom_bed_file
# -------------------------------------------------------

function main_func {
	i=$1
	zcat < $cov_file_dir/gnomad.genomes.r${version}.chr${i}.coverage.txt.gz | grep -v '#' | awk -v d="$min_depth" '{if($3 > d) {print $1"\t"$2-1"\t"$2"\t"$3}}' > $out_dir/gnomad.chr${i}.high_cov.single_nt.bed
	mergeBed -i $out_dir/gnomad.chr${i}.high_cov.single_nt.bed > $out_dir/gnomad.chr${i}.high_cov.bed
}


# ============== Main Run ==============
for i in `seq 1 22` X Y
do
	echo $i
	main_func $i &
done

wait
# ======================================




# -------- Merge bed files from each chromosome into a single BED file --------
for i in `seq 1 22` X Y
do
	cat $out_dir/gnomad.chr${i}.high_cov.bed >> $out_all_chrom_bed_file
done

# -----------------------------------------------------------------------------

# ---- Calculate total coverage ----
cat $out_all_chrom_bed_file | awk '{sum+=$3-$2} END {print sum}' > $out_dir/Total_coverage.min_depth${min_depth}.txt
