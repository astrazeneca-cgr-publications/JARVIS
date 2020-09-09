#!/bin/bash
#SBATCH -o out.subset_scores
#SBATCH -e err.subset_scores
#SBATCH -t 24:0:0
#SBATCH -n 4
#SBATCH --mem=8G


win_len=$1

sv_bed_dir="sv-bed"
sv_classes=('utr' 'intergenic' 'lof' 'promoter' 'copy_gain' 'dup_partial' 'inv_span' 'dup_lof' 'intronic')



# ===================  phastCons primate  ===================
subset_phastcons_primate_across_all_chr () {

	cl=$1
	echo ">> SV class: $cl"
	
	phastcons_base_dir="/projects/cgr/users/kclc950/JARVIS/other_datasets/conservation/phastCons46way_primates/bed"
	out_scores_dir="phastCons46way-sv-bed"
	mkdir -p $out_scores_dir

	for chr in `seq 1 22`; do
		echo "Chr: $chr"
		chr_phastcons_file="$phastcons_base_dir/chr${chr}.phastCons46way.primates.high_conf_regions.bed.gz"

		tabix $chr_phastcons_file -B $sv_bed_dir/SV.${cl}.bed > $out_scores_dir/phastCons46way.SV.${cl}.chr${chr}.bed &
	done

	wait

	cat $out_scores_dir/phastCons46way.SV.${cl}.chr*.bed > $out_scores_dir/phastCons46way.SV.${cl}.with_coords.bed
	rm $out_scores_dir/phastCons46way.SV.${cl}.chr*.bed

	echo ">> SV class: $cl -- DONE"
}



subset_phastcons_primate () {
	echo -e "\n> phastCons-primate"

	for cl in "${sv_classes[@]}"; do
		subset_phastcons_primate_across_all_chr $cl &
	done
	wait
}








# ===================  gwRVIS  ===================
subset_gwrvis_across_all_chr () {

	cl=$1
	echo ">> SV class: $cl"
	
	gwrvis_base_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-single-nt-gwrvis-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/single_nt_gwrvis_per_chr/"
	out_scores_dir="gwrvis-sv-bed"
	mkdir -p $out_scores_dir

	for chr in `seq 1 22`; do
	#for chr in `seq 21 21`; do
		echo "Chr: $chr"
		chr_gwrvis_file="$gwrvis_base_dir/gwrvis_single_nt.chr${chr}.bed.gz"

		tabix $chr_gwrvis_file -B $sv_bed_dir/SV.${cl}.bed > $out_scores_dir/gwrvis.SV.${cl}.chr${chr}.bed &
	done
	wait

	cat $out_scores_dir/gwrvis.SV.${cl}.chr*.bed > $out_scores_dir/gwrvis.SV.${cl}.with_coords.bed
	rm $out_scores_dir/gwrvis.SV.${cl}.chr*.bed

	echo ">> SV class: $cl -- DONE"
}


subset_gwrvis () {
	echo -e "\n> gwrvis"

	for cl in "${sv_classes[@]}"; do
		subset_gwrvis_across_all_chr $cl &
	done
	wait
}












# ===================  JARVIS  ===================
subset_jarvis_across_all_chr () {

	cl=$1
	echo ">> SV class: $cl"
	
	jarvis_base_dir="/projects/cgr/users/kclc950/JARVIS/out/topmed-NEW_genome_wide_scores-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/jarvis_predictions/"
	out_scores_dir="jarvis-sv-bed"
	mkdir -p $out_scores_dir

	for chr in `seq 1 22`; do
	#for chr in `seq 21 21`; do
		echo "Chr: $chr"
		chr_jarvis_file="$jarvis_base_dir/chr${chr}/jarvis.chr${chr}.both-features.sorted.bed.bgz"

		tabix $chr_jarvis_file -B $sv_bed_dir/SV.${cl}.${no_chr_str}bed > $out_scores_dir/jarvis.SV.${cl}.chr${chr}.bed &
	done
	wait

	cat $out_scores_dir/jarvis.SV.${cl}.chr*.bed > $out_scores_dir/jarvis.SV.${cl}.with_coords.bed
	rm $out_scores_dir/jarvis.SV.${cl}.chr*.bed

	echo ">> SV class: $cl -- DONE"
}


subset_jarvis () {
	echo -e "\n> jarvis"

	for cl in "${sv_classes[@]}"; do
		subset_jarvis_across_all_chr $cl &
	done
	wait
}




add_af_to_jarvis () {

	out_scores_dir="jarvis-sv-bed"

	for cl in "${sv_classes[@]}"; do
		echo ">> $cl"	
		intersectBed -wo -a $out_scores_dir/jarvis.SV.${cl}.with_coords.bed -b sv-bed/SV.${cl}.bed | awk '{print $8"\t"$4}' > $out_scores_dir/jarvis.SV.${cl}.bed &
	done
	wait
}







intersect_wrapper () {
	score_ref_file=$1
	cl=$2
	sv_bed_dir=$3
	out_scores_dir=$4
	score=$5	
	
	no_chr_str=""
	if [ $score == "cadd" ] || [ $score == "orion" ]; then
		no_chr_str="no_chr."
	fi	


	tabix $score_ref_file -B $sv_bed_dir/SV.${cl}.${no_chr_str}bed > $out_scores_dir/${score}.SV.${cl}.bed.tmp
	return 0  # TEMP


	if [ $score == "linsight" ] || [ $score == "ncER_10bp" ]; then
		mv $out_scores_dir/${score}.SV.${cl}.bed.tmp $out_scores_dir/${score}.SV.${cl}.with_coords.bed

		intersectBed -wo -a $out_scores_dir/${score}.SV.${cl}.bed.tmp -b $sv_bed_dir/SV.${cl}.${no_chr_str}bed | awk '{print $8"\t"$4}' > $out_scores_dir/${score}.SV.${cl}.bed

	elif [ $score == "cadd" ]; then
		cat $out_scores_dir/${score}.SV.${cl}.bed.tmp | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > $out_scores_dir/${score}.SV.${cl}.bed.tmp2; mv $out_scores_dir/${score}.SV.${cl}.bed.tmp2 $out_scores_dir/${score}.SV.${cl}.bed.tmp

		intersectBed -wo -a $out_scores_dir/${score}.SV.${cl}.bed.tmp -b $sv_bed_dir/SV.${cl}.${no_chr_str}bed | awk '{print $8"\t"$4}' > $out_scores_dir/${score}.SV.${cl}.bed

	elif [ $score == "orion" ]; then
		intersectBed -wo -a $out_scores_dir/${score}.SV.${cl}.bed.tmp -b $sv_bed_dir/SV.${cl}.${no_chr_str}bed | awk '{print $9"\t"$4}' > $out_scores_dir/${score}.SV.${cl}.bed
	fi

	#rm $out_scores_dir/${score}.SV.${cl}.bed.tmp
}



subset_linsight () {
	echo -e "\n> linsight"
	score_ref_file="../../other_datasets/genome-wide-scores/linsight/linsight.bed.bgz"
	out_scores_dir="linsight-sv-bed"
	mkdir -p $out_scores_dir

	for cl in "${sv_classes[@]}"; do
		echo $cl
		intersect_wrapper $score_ref_file $cl $sv_bed_dir $out_scores_dir linsight &
	done
}


subset_ncER_10bp () {
	echo -e "\n> ncER_10bp"
	score_ref_file="../../other_datasets/genome-wide-scores/ncER_10bp/ncER_10bpBins_allChr_coordSorted.txt.bgz"
	out_scores_dir="ncER_10bp-sv-bed"
	mkdir -p $out_scores_dir

	for cl in "${sv_classes[@]}"; do
		echo $cl
		intersect_wrapper $score_ref_file $cl $sv_bed_dir $out_scores_dir ncER_10bp &
	done
}


subset_cadd () {
	echo -e "\n> cadd"
	score_ref_file="../../other_datasets/genome-wide-scores/cadd/cadd.whole_genome.all_raw_scores.vcf.bgz"
	out_scores_dir="cadd-sv-bed"
	mkdir -p $out_scores_dir

	for cl in "${sv_classes[@]}"; do
		intersect_wrapper $score_ref_file $cl $sv_bed_dir $out_scores_dir cadd &
	done
}


subset_orion () {
	echo -e "\n>orion"
	score_ref_file="../../other_datasets/genome-wide-scores/orion/orion.1001.masked.new.txt.gz"
	out_scores_dir="orion-sv-bed"
	mkdir -p $out_scores_dir

	for cl in "${sv_classes[@]}"; do
		intersect_wrapper $score_ref_file $cl $sv_bed_dir $out_scores_dir orion &
	done
}


#subset_jarvis
#add_af_to_jarvis
subset_gwrvis

#subset_phastcons_primate
#subset_linsight
#subset_ncER_10bp
#subset_cadd
#subset_orion

wait
