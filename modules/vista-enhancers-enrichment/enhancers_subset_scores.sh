#!/bin/bash
#SBATCH -o out.subset_scores
#SBATCH -e err.subset_scores
#SBATCH -t 24:0:0
#SBATCH -n 4
#SBATCH --mem=4G


input_bed=$1

# TODO: Replace these with your own PATHs
jarvis_base_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-NEW_genome_wide_scores-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/jarvis_predictions/" 
gwrvis_base_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-single-nt-gwrvis-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/single_nt_gwrvis_per_chr/" 




# ===================  gwRVIS  ===================
subset_gwrvis () {

	input_bed=$1
	echo ">> gwRVIS - Input: $input_bed"
	

	for chr in `seq 1 22`; do
		echo "Chr: $chr"
		chr_gwrvis_file="$gwrvis_base_dir/gwrvis_single_nt.chr${chr}.bed.gz"

		tabix $chr_gwrvis_file -B $input_bed > ${input_bed}.gwrvis.chr${chr}.tmp &
	done
	wait

	cat ${input_bed}.gwrvis.chr*.tmp > ${input_bed}.gwrvis.bed 
	rm ${input_bed}.gwrvis.chr*.tmp

	echo ">> gwRVIS - Input: $input_bed -- [DONE]"
}





# ===================  JARVIS  ===================
subset_jarvis () {

	input_bed=$1
	echo ">> JARVIS - Input: $input_bed"
	

	for chr in `seq 1 22`; do
		echo "Chr: $chr"
		chr_jarvis_file="$jarvis_base_dir/chr${chr}/jarvis.chr${chr}.both-features.sorted.bed.bgz"

		tabix $chr_jarvis_file -B $input_bed > ${input_bed}.jarvis.chr${chr}.tmp &
	done
	wait

	cat ${input_bed}.jarvis.chr*.tmp > ${input_bed}.jarvis.bed 
	rm ${input_bed}.jarvis.chr*.tmp

	echo ">> JARVIS - Input: $input_bed -- [DONE]"
}





# ========== MAIN RUN ==========
subset_jarvis $input_bed &
subset_gwrvis $input_bed &
wait
