#!/bin/bash
#SBATCH -t 24:0:0
#SBATCH -n 1
#SBATCH --mem=8G

sv_bed_dir="sv-bed"
sv_classes=('intergenic' 'lof' 'promoter' 'copy_gain' 'utr' 'dup_partial' 'inv_span' 'dup_lof' 'intronic')


subset_jarvis () {
	echo -e "\n> jarvis"
	jarvis_ref_file="../../jarvis-gwrvis-scores/bed/JARVIS.prediction_scores.bed"
	out_scores_dir="jarvis-sv-bed"
	mkdir -p $out_scores_dir

	for cl in "${sv_classes[@]}"; do
		echo $cl
		intersectBed -wo -a $sv_bed_dir/SV.${cl}.bed -b $jarvis_ref_file | cut -f4,9,10 > $out_scores_dir/jarvis.SV.${cl}.bed
	done
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


	if [ $score == "linsight" ] || [ $score == "ncER_10bp" ]; then
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
#subset_linsight
subset_ncER_10bp
#subset_cadd
#subset_orion

wait
