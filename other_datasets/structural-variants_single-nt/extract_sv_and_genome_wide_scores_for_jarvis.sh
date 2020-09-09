#!/bin/bash
#SBATCH -o out.compare_sv_with_genome_wide_distr.jarvis
#SBATCH -n 4
#SBATCH --time=12:00:00


# ==== Input dirs & files ====
# JARVIS genome wide scores
jarvis_base_dir="../../out/topmed-NEW_genome_wide_scores-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/jarvis_predictions"

# - genome wide BED
classes_bed_dir="../../out/topmed-NEW_genome_wide_scores-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/gwrvis_distribution"

# SV BED
sv_bed_dir="../../other_datasets/structural-variants_single-nt/sv-bed"


# ==== Output dirs & files ====
out_dir='gen-pop_vs_SV-distr'
mkdir -p $out_dir







# Get JARVIS scores for UTRs
function subset_score_for_chr () {

	ref_bed_file=$1
	class=$2 #utr or intronic/intron
	chr=$3
	file_identifier=$4

	echo "Class: $class - Chr $chr"
	cur_jarvis_bed="$jarvis_base_dir/chr$chr/jarvis.chr${chr}.both-features.sorted.bed.bgz"

	tabix $cur_jarvis_bed -B $ref_bed_file | cut -f4 > $out_dir/jarvis.${file_identifier}${class}.chr${chr}.tmp
	echo "Class: $class - Chr $chr -- [DONE]"
}



function process_class_for_jarvis () {

	ref_bed_file=$1
	class=$2
	file_identifier=$3

	echo "ref_bed_file: $ref_bed_file"

	for chr in `seq 1 22`; do
	#for chr in `seq 21 21`; do
		subset_score_for_chr $ref_bed_file $class $chr $file_identifier &
	done
	wait

	cat $out_dir/jarvis.${file_identifier}${class}.chr*.tmp > $out_dir/jarvis.${file_identifier}${class}.txt
	rm $out_dir/jarvis.${file_identifier}${class}.chr*.tmp
}




# ============ RUN for JARVIS ============
# UTR
class=utr
file_identifier="sv_"
sv_class_bed="$sv_bed_dir/SV.${class}.bed"
process_class_for_jarvis $sv_class_bed $class $file_identifier &

file_identifier=""
class_bed="$classes_bed_dir/gwrvis.${class}.mutually_excl.bed"
process_class_for_jarvis $class_bed $class $file_identifier &



# Intron
#class=intronic
#file_identifier="sv_"
#sv_class_bed="$sv_bed_dir/SV.${class}.bed"
#process_class_for_jarvis $sv_class_bed $class $file_identifier &

#class=intron
#file_identifier=""
#class_bed="$classes_bed_dir/gwrvis.${class}.mutually_excl.bed"
#process_class_for_jarvis $class_bed $class $file_identifier &
wait
