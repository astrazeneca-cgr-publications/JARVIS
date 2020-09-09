#!/bin/bash
#SBATCH -o out.compare_sv_with_genome_wide_distr.gwrvis
#SBATCH -n 4
#SBATCH --time=12:00:00


# ==== Input dirs & files ====
# gwRVIS genome wide scores
gwrvis_base_dir="../../out/topmed-single-nt-gwrvis-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/single_nt_gwrvis_per_chr"

# - genome wide BED
classes_bed_dir="../../out/topmed-NEW_genome_wide_scores-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/gwrvis_distribution"

# SV BED
sv_bed_dir="../../other_datasets/structural-variants_single-nt/sv-bed"


# ==== Output dirs & files ====
out_dir='gen-pop_vs_SV-distr'
mkdir -p $out_dir





# Get gwRVIS scores for UTRs
function subset_score_for_chr () {

	ref_bed_file=$1
	class=$2 #utr or intronic/intron
	chr=$3
	file_identifier=$4

	echo "Class: $class - Chr $chr"
	cur_gwrvis_bed="$gwrvis_base_dir/gwrvis_single_nt.chr${chr}.bed.gz"

	tabix $cur_gwrvis_bed -B $ref_bed_file | cut -f4 > $out_dir/gwrvis.${file_identifier}${class}.chr${chr}.tmp
	echo "Class: $class - Chr $chr -- [DONE]"
}



function process_class_for_gwrvis () {

	ref_bed_file=$1
	class=$2
	file_identifier=$3

	echo "ref_bed_file: $ref_bed_file"

	for chr in `seq 1 22`; do
	#for chr in `seq 21 21`; do
		subset_score_for_chr $ref_bed_file $class $chr $file_identifier &
	done
	wait

	cat $out_dir/gwrvis.${file_identifier}${class}.chr*.tmp > $out_dir/gwrvis.${file_identifier}${class}.txt
	rm $out_dir/gwrvis.${file_identifier}${class}.chr*.tmp
}







# ============ RUN for gwRVIS ============
# UTR
class=utr
file_identifier="sv_"
sv_class_bed="$sv_bed_dir/SV.${class}.bed"
process_class_for_gwrvis $sv_class_bed $class $file_identifier &

file_identifier=""
class_bed="$classes_bed_dir/gwrvis.${class}.mutually_excl.bed"
process_class_for_gwrvis $class_bed $class $file_identifier &



# Intron
class=intronic
file_identifier="sv_"
sv_class_bed="$sv_bed_dir/SV.${class}.bed"
process_class_for_gwrvis $sv_class_bed $class $file_identifier &

class=intron
file_identifier=""
class_bed="$classes_bed_dir/gwrvis.${class}.mutually_excl.bed"
process_class_for_gwrvis $class_bed $class $file_identifier &
wait
