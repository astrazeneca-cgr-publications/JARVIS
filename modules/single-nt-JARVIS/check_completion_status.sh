#!/bin/bash

out_dir="../../out/topmed-genome_wide_scores_v2-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/jarvis_predictions/"


function check_chr {

	chr=$1
	chr_dir="$out_dir/chr$chr"

	if [ -d "$chr_dir" ]; then
		echo -e "\nProcessing Chr $chr..."

		chr_interv_files=`ls $chr_dir/chr${chr}.both-features.* | wc -l`
		covered_real_estate=$(($chr_interv_files*50000))

		jarvis_pred=`ls $chr_dir/chr${chr}.both-features.* | sed 's/.*\///g' |sed "s/chr${chr}.both-features.//g" | sed 's/.jarvis//g' | sort -V | tail -1`
		jarvis_pred_right=`echo $jarvis_pred | sed "s/.*\_//"`


		real_est=`zcat < ../../single-nt-jarvis-gwrvis-scores/JARVIS/jarvis.chr${chr}.both-features.sorted.bed.bgz | wc -l`
		
		status='Pending...'
		if [ $covered_real_estate -eq $jarvis_pred_right -a $jarvis_pred_right -ge $real_est ]; then
			status="[DONE] - chr${chr}"
		fi


		full_str="\n> Chr: $chr"
		full_str="${full_str}\nJarvis predictions: $jarvis_pred"
		full_str="${full_str}\nFull real estate: $real_est"
		full_str="${full_str}\nCovered real estate (est.): $covered_real_estate"
		full_str="${full_str}\n--Status: $status"
		echo -e $full_str
	fi
}


for chr in `seq 1 22`; do
#for chr in 21; do 
	check_chr $chr & 	
done
wait
