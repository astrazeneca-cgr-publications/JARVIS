#!/bin/bash 
#SBATCH -o out.run_wgs_cov_qc.txt
#SBATCH -e err.run_wgs_cov_qc.txt
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G 
#SBATCH --time=24:00:00

out_dir="../coverage-files"
out_file="$out_dir/All_mean_depths.gnomad.WGS_coverage.txt"

function concat_all_files {
	rm -f $out_file

	for i in `seq 1 22`
	do
		zcat < $out_dir/gnomad.genomes.r2.0.2.chr${i}.coverage.txt.gz | grep -v '#' | awk '{print $3}' >> $out_file
	done

	i='X'
	zcat < $out_dir/gnomad.genomes.r2.0.2.chr${i}.coverage.txt.gz | grep -v '#' | awk '{print $3}' >> $out_file
}

# > Concatenate mean depths per nt from all chromosomes
concat_all_files


# > Sort file with mean depths per nt across all chromosomes -- Reuqired only to extract percentile metrics, e.g. median
sort -n --parallel=8 $out_file > ${out_file}.sorted
# > Calculate median depth across all chromosomes
median=`cat ${out_file}.sorted | awk '{vals[NR] = $1} END{if (NR % 2) {print vals[(NR + 1) / 2]} else { print(vals[(NR / 2)] + vals[(NR / 2) + 1]) / 2.0}}'`


# > Calculate min, max and mean depth across all chromosomes and save all metrics to a file
cat $out_file | awk -v median="$median" 'NR == 1 || $1 < min {min = $1}; NR == 1 || $1 > max {max = $1}; {sum+=$1} END {print "Min: " min; print "Max: " max; print "Mean: " sum / NR; print "Median: " median}' > $out_dir/WGS_coverage.QC_depths.txt
