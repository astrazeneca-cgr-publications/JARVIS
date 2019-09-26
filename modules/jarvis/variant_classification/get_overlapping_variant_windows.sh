#!/bin/sh

# TODO: read these from input_classes.txt file
genomic_classes_file=$1 
out_dir=$2 
clinvar_feature_table_dir=$3 
patho_benign_sets=$4

#genomic_classes_file=/projects/cgr/users/kclc950/gwRVIS/modules/genomic_classes.log 
#out_dir=/projects/cgr/users/kclc950/gwRVIS/out/gnomad-regression_beta-winlen_3000.MAF_0.001.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions 
#clinvar_feature_table_dir=/projects/cgr/users/kclc950/gwRVIS/out/gnomad-regression_beta-winlen_3000.MAF_0.001.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions/ml_data/clinvar_feature_tables


declare -a genomic_classes=(`tail -n+2 $genomic_classes_file | cut -f1`)
declare -a genomic_classes_files=(`tail -n+2 $genomic_classes_file | cut -f2`)


declare -A input_classes
for i in `seq 1 $((${#genomic_classes[@]}-1))`; do   
	#echo "${genomic_classes[$i]}" 
	input_classes["${genomic_classes[$i]}"]="${genomic_classes_files[$i]}"
done
#for i in "${!input_classes[@]}"; do
#	printf "$i\t${input_classes[$i]}\n"
#done

coverage_ratio_diff_thres=0.15

# -----------------------------------------------------

# > Create BED with all gwRVIS windows
tail -n+2 $out_dir/gwrvis_scores/full_genome.all_gwrvis.bed | cut -f2,3,4 > $clinvar_feature_table_dir/all_gwrvis_windows.bed


# > Get full feature tables per genomic class
for gen_class in ccds utr lincrna ucne intergenic; do
	tail -n+2 $clinvar_feature_table_dir/full_feature_table.${patho_benign_sets}.bed | awk -v gen_class="$gen_class" '{if($6 == gen_class) print $0}' > $clinvar_feature_table_dir/full_feature_table.${patho_benign_sets}.${gen_class}.bed
done


# > Create BED files: full gwRVIS windows for CCDS only; region-specific BED files for all other non-coding regions
intersectBed -wa -a $clinvar_feature_table_dir/all_gwrvis_windows.bed -b $clinvar_feature_table_dir/full_feature_table.${patho_benign_sets}.ccds.bed | sort -k1,1 -k2,2n | mergeBed > $clinvar_feature_table_dir/ccds_gwrvis_windows.bed


for gen_class in ccds utr lincrna ucne intergenic; do
	cat $clinvar_feature_table_dir/full_feature_table.${patho_benign_sets}.${gen_class}.bed | cut -f1,2,3 > $clinvar_feature_table_dir/${gen_class}_gwrvis_win_subregions.bed
done



for overlap_class in intergenic utr; do
	
	echo $overlap_class

	# Get BED file with all gwRVIS windows that CCDS and the 'overlap_class' have variants, including pathogenic/benign content annotation
	intersectBed -wa -a $clinvar_feature_table_dir/ccds_gwrvis_windows.bed -b $clinvar_feature_table_dir/${overlap_class}_gwrvis_win_subregions.bed |  intersectBed -wao -a stdin -b $clinvar_feature_table_dir/full_feature_table.${patho_benign_sets}.${overlap_class}.bed | awk '{print $1"\t"$2"\t"$3"\t"$7}' | mergeBed -c 4 -o distinct > $clinvar_feature_table_dir/${overlap_class}.overlap_with_ccds.bed

	# Get real-estate of overlapping windows covered by CCDS
	intersectBed -a $clinvar_feature_table_dir/${overlap_class}.overlap_with_ccds.bed -b "${input_classes['ccds']}" | awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' | intersectBed -wo -b stdin -a $clinvar_feature_table_dir/ccds_gwrvis_windows.bed | cut -f1,2,3,7,8 | mergeBed -c 4,5 -o sum,distinct > $clinvar_feature_table_dir/ccds_coverage.${overlap_class}_overlaping_ccds.bed


	# Get real-estate of overlapping windows covered by the 'overlap_class'
	intersectBed -a $clinvar_feature_table_dir/${overlap_class}.overlap_with_ccds.bed -b "${input_classes[$overlap_class]}"  | awk '{print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4}' | intersectBed -wo -b stdin -a $clinvar_feature_table_dir/ccds_gwrvis_windows.bed | cut -f1,2,3,7,8 | mergeBed -c 4,5 -o sum,distinct > $clinvar_feature_table_dir/${overlap_class}_coverage.${overlap_class}_overlaping_ccds.bed


	intersectBed -wo -a $clinvar_feature_table_dir/ccds_coverage.${overlap_class}_overlaping_ccds.bed -b $clinvar_feature_table_dir/${overlap_class}_coverage.${overlap_class}_overlaping_ccds.bed | awk '{print $1"\t"$2"\t"$3"\t"$4/3000"\t"$9/3000"\t"($4-$9)/3000"\t"$10 }' |  awk -v diff_thres="$coverage_ratio_diff_thres" '{ if($6 > diff_thres) print $0}' > $clinvar_feature_table_dir/blacklist.${overlap_class}_overlaping_ccds.bed
	# Output columns are as follows: chr     start   end     ccds_cov        other_class_cov  diff (CCDS-other_class)    variant_type (pathogenic/1 or benign/0)

done
