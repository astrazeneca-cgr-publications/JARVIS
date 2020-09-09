#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --mem=16GB

##SBATCH --cpus-per-task=2
#SBATCH --partition=gpu
#SBATCH --gres=gpu:volta:1


win_len=$1

# Global variables
ccds_bed_base_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-CCDS-ClinVar_pathogenic-denovodb_benign-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/gwrvis_distribution"
declare -a ccds_classes=("omim-HI" "intolerant" "tolerant" "ccds")


single_nt_gwrvis_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-single-nt-gwrvis-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/single_nt_gwrvis_per_chr"

out_base_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-single-nt-gwrvis-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/gwrvis_distr_in_ccds_classes"
mkdir -p $out_base_dir


function subset_gwrvis_for_ccds_class {

	cl=$1

	cur_ccds_bed="${ccds_bed_base_dir}/gwrvis.${cl}.mutually_excl.bed"
	echo "CCDS class: $cl"

	cur_ccds_out_dir="$out_base_dir/$cl"
	mkdir -p $cur_ccds_out_dir


	for chrom in `seq 1 22`; do
		echo "chr: $chrom"
		cur_gwrvis_bed="${single_nt_gwrvis_dir}/gwrvis_single_nt.chr${chrom}.bed.gz"

		tabix $cur_gwrvis_bed -B $cur_ccds_bed > $cur_ccds_out_dir/temp.chr${chrom}.bed &
	done
	wait


	cat $cur_ccds_out_dir/temp.chr*.bed > $cur_ccds_out_dir/${cl}.gwrvis.bed
}



# === Main run ===
for cl in "${ccds_classes[@]}"; do
	subset_gwrvis_for_ccds_class $cl &
done

wait

