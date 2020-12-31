#!/bin/bash
#SBATCH -o out.compile_single_nt_jarvis_per_chr
#SBATCH -e err.compile_single_nt_jarvis_per_chr
#SBATCH -n 8
#SBATCH --mem=120G
#SBATCH --time=24:00:00

win_len=3000
base_jarvis_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-genome_wide_scores_v2-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/jarvis_predictions"


declare -a complete_chroms=(16 17 18 20 21 22)


function concat_jarvis_files {
	echo "Concatenating jarvis files per chrom..."
	
	for chrom in "${complete_chroms[@]}"; do
		echo "Chr: $chrom"
		jarvis_chrom_dir="${base_jarvis_dir}/chr${chrom}"

		cat $jarvis_chrom_dir/chr${chrom}.both-features.*.jarvis > $jarvis_chrom_dir/jarvis.chr${chrom}.both-features.bed &
	done
	wait
}



function sort_jarvis_files {
	echo "Sorting jarvis files per chrom..."

	cnt=1
	for chrom in "${complete_chroms[@]}"; do
		echo "Chr: $chrom"
		jarvis_chrom_dir="${base_jarvis_dir}/chr${chrom}"

		sort -V -k 1,1 -k2,2n $jarvis_chrom_dir/jarvis.chr${chrom}.both-features.bed > $jarvis_chrom_dir/jarvis.chr${chrom}.both-features.sorted.bed &

		cnt=$((cnt+1))
		if [ "$cnt" == 10 ]; then
			cnt=1
			wait
		fi
	done
	wait
}



function bgzip_sorted_jarvis_files {
	echo "Bgzip sorted jarvis files per chrom..."

	cnt=1
	for chrom in "${complete_chroms[@]}"; do
		echo "Chr: $chrom"
		jarvis_chrom_dir="${base_jarvis_dir}/chr${chrom}"

		cat $jarvis_chrom_dir/jarvis.chr${chrom}.both-features.sorted.bed | bgzip > $jarvis_chrom_dir/jarvis.chr${chrom}.both-features.sorted.bed.bgz &

		#cnt=$((cnt+1))
		#if [ "$cnt" == 10 ]; then
		#	cnt=1
		#	wait
		#fi
	done
	wait
}




function build_tabix_for_jarvis_files {
	echo "Build tabix for jarvis files per chrom..."

	cnt=1
	for chrom in "${complete_chroms[@]}"; do
		echo "Chr: $chrom"
		jarvis_chrom_dir="${base_jarvis_dir}/chr${chrom}"

		tabix -p bed $jarvis_chrom_dir/jarvis.chr${chrom}.both-features.sorted.bed.bgz &

		#cnt=$((cnt+1))
		#if [ "$cnt" == 10 ]; then
		#	cnt=1
		#	wait
		#fi
	done
	wait
}



# ==== MAIN ====
concat_jarvis_files


sort_jarvis_files


bgzip_sorted_jarvis_files


build_tabix_for_jarvis_files
