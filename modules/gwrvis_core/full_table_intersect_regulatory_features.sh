#!/bin/bash

# Append regulatory features with single-nt resolution to the full feature table BED file
full_feature_table=$1 #"full_feature_table.clinvar_denovodb.bed"



# First delete previous window-based annotation
cat $full_feature_table | cut -f1-19 > ${full_feature_table}.tmp
mv ${full_feature_table}.tmp $full_feature_table




cell_line="Monocytes_CD14plus"

for regul_elem in CTCF_Binding_Site Enhancer Open_chromatin TF_binding_site H3K27ac H3K27me3 H4K20me1 H3K9ac H3K4me1 H3K4me2 H3K4me3 H3K36me3; do
#for regul_elem in Enhancer CTCF_Binding_Site Open_chromatin TF_binding_site H3K27ac H3K27me3 H4K20me1 H3K9ac H3K4me1 H3K4me2 H3K4me3 H3K36me3; do

	echo -e "\n\n>>Regulatory element: $regul_elem"

	ensembl_regul_file="../ensembl/GRCh37-Regulatory_Features/${cell_line}.${regul_elem}.sorted.bed"
	echo $ensembl_regul_file
	
	tmp_regul_inters_bed=${full_feature_table}_${regul_elem}.tmp
	# Add column name for current regulatory element
	echo $regul_elem > $tmp_regul_inters_bed

	
	# record 1/0 based on intersection or not
	tail -n+2 $full_feature_table | intersectBed -wao -a stdin -b $ensembl_regul_file |  awk '!seen[$1"_"$2] {print} {++seen[$1"_"$2]}' | rev | cut -f1 | rev | awk '{if($1 > 0) {print 1} else print 0}' >> $tmp_regul_inters_bed
	cat $tmp_regul_inters_bed | wc -l

	
	# append 1/0 to the full feature file
	paste $full_feature_table $tmp_regul_inters_bed > ${full_feature_table}.tmp
	mv ${full_feature_table}.tmp ${full_feature_table}

	echo ${full_feature_table}


	# clean up
	rm $tmp_regul_inters_bed
	
done
