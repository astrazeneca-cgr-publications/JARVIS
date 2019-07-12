#!/bin/sh

out_dir=$1
input_classes_file=$2
readarray -t input_classes < $input_classes_file

gwrvis_distr_dir="${out_dir}/gwrvis_distribution"
full_feature_table="${out_dir}/full_genome_out/full_gwrvis_and_regulatory_features.bed"

ml_data_dir="${out_dir}/ml_data"
mkdir -p  $ml_data_dir

out_feature_table_dir=$ml_data_dir/feature_tables
mkdir -p $out_feature_table_dir

# Output fle
out_full_feature_table="${out_feature_table_dir}/full_gwrvis_and_regulatory_features.All_genomic_classes.tsv"


function add_header_to_bed_by_genomic_class () {

	full_feature_table_by_genomic_class=$1

	full_header=()
	tmp_header=( $(cat $full_feature_table | head -1) )

	# add 'genomic_class' column right after 'gwrvis' column
	col_after_gwrvis=false
	for v in "${tmp_header[@]}"; do
		if [ $col_after_gwrvis = true ]; then
			full_header+=( "genomic_class" )
			col_after_gwrvis=false
		fi
		full_header+=( "$v" )

		if [ $v == 'gwrvis' ]; then
			col_after_gwrvis=true
		fi
	done

	printf "%s\t" "${full_header[@]}" > ${full_feature_table_by_genomic_class}.header
	sed -i 's/$/\n/' ${full_feature_table_by_genomic_class}.header
	cat ${full_feature_table_by_genomic_class}.header ${full_feature_table_by_genomic_class}.tmp > ${full_feature_table_by_genomic_class}
	rm ${full_feature_table_by_genomic_class}.header ${full_feature_table_by_genomic_class}.tmp
}




function get_feature_table_by_genomic_class {

	for class in "${input_classes[@]}"; do
		echo $cl
		genomic_class_out_file="${gwrvis_distr_dir}/gwrvis.${class}.mutually_excl.bed"

		echo "Getting full genomic class BED file across all chromosomes ... (mutually exclusive with other classes)"
		cat ${gwrvis_distr_dir}/BED/gwrvis_scores_chr*.genomic_coords.${class}.bed | sed "s/\$/\t$class/g" > $genomic_class_out_file
		echo "[$class] Total nt length: "`cat $genomic_class_out_file | awk '{sum+=$3-$2} END {print sum}'`


		echo -e "Retrieving feature table for $class ...\n"
		full_feature_table_by_genomic_class=${out_feature_table_dir}/full_gwrvis_and_regulatory_features.${class}.tsv



		# Keep coordinates from the genomic class regions (not the full 3kb windows) -- remove regions with 'NaN' gwRVIS value
		tail -n +2 $full_feature_table | intersectBed -wao -a $genomic_class_out_file -b stdin| grep $class | grep -v "NaN" | sortBed | cut --complement -f6,7,8,9,37 > ${full_feature_table_by_genomic_class}.tmp


		# adding header to feature table by class and saving into file
		add_header_to_bed_by_genomic_class $full_feature_table_by_genomic_class
	done
}


function merge_feature_tables_from_all_classes {

	echo "Merging feature tables from all genomic classes ..."
	rm -f $out_full_feature_table 

	# Minor bug: remove a single overlapping small region (between an intron and a lincrna) -- calling awk !seen for that ...
	cat ${out_feature_table_dir}/full_gwrvis_and_regulatory_features.*.tsv | grep -v gwrvis | awk '!seen[$1"_"$2] {print} {++seen[$1"_"$2]}' > ${out_full_feature_table}.tmp

	# add header
	cat ${out_feature_table_dir}/full_gwrvis_and_regulatory_features.*.tsv | head -1 > ${out_full_feature_table}.header
	cat ${out_full_feature_table}.header ${out_full_feature_table}.tmp > $out_full_feature_table
	rm ${out_full_feature_table}.tmp ${out_full_feature_table}.header
}



function add_clinvar_annotation {

	clinvar_feature_table_dir=$ml_data_dir/clinvar_feature_tables
	mkdir -p $clinvar_feature_table_dir

	# get intersections with clinvar pathogenic/benign
	tail -n +2 $out_full_feature_table | intersectBed -wo -a ../other_datasets/clinvar/clinvar.pathogenic.bed -b stdin | cut --complement -f5,6,7,37 > $clinvar_feature_table_dir/full_feature_table.clinvar_pathogenic.bed.tmp
	tail -n +2 $out_full_feature_table | intersectBed -wo -a ../other_datasets/clinvar/clinvar.benign.bed -b stdin | cut --complement -f5,6,7,37 > $clinvar_feature_table_dir/full_feature_table.clinvar_benign.bed.tmp

	cat $clinvar_feature_table_dir/full_feature_table.clinvar_pathogenic.bed.tmp $clinvar_feature_table_dir/full_feature_table.clinvar_benign.bed.tmp > $clinvar_feature_table_dir/full_feature_table.clinvar.bed.tmp

	# add header
	tmp_header_file=$clinvar_feature_table_dir/header.tmp
	cat $out_full_feature_table | head -1 | sed 's/gwrvis/clinvar_annot\tgwrvis/' | sed 's/\t$//' > $tmp_header_file
	cat $tmp_header_file $clinvar_feature_table_dir/full_feature_table.clinvar.bed.tmp > $clinvar_feature_table_dir/full_feature_table.clinvar.bed

	# clenaup tmp files
	#rm $clinvar_feature_table_dir/*.tmp
}


# =============== MAIN RUN ================

# Merge BED files by genomic class and retrieve respective feature table
#get_feature_table_by_genomic_class

# merge feature tables from all genomic classes
#merge_feature_tables_from_all_classes


# annotate full feature table (already with genomic class annotation) with ClinVar pathogenic/bengign variants
add_clinvar_annotation
