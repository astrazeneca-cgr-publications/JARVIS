#!/bin/sh

# ==== Read / Infer input arguments ====
config_file=$1
input_classes_file=$2
readarray -t input_classes < $input_classes_file

out_dir=`python custom_utils.py $config_file`
printf "Out dir: ${out_dir}\n"


pathogenic_set=`cat $config_file | grep "pathogenic_set:" | sed 's/^[ \t]*pathogenic_set: //' | sed 's/[ ]*$//'`
benign_set=`cat $config_file | grep "benign_set:" | sed 's/^[ \t]*benign_set: //' | sed 's/[ ]*$//'`
echo "Pathogenic set: $pathogenic_set"
echo "Benign set: $benign_set"



# ==== Global variables ====
gwrvis_distr_dir="${out_dir}/gwrvis_distribution"
full_feature_table="${out_dir}/full_genome_out/full_gwrvis_and_regulatory_features.bed"

ml_data_dir="${out_dir}/ml_data"
mkdir -p  $ml_data_dir

out_feature_table_dir=$ml_data_dir/feature_tables
mkdir -p $out_feature_table_dir

clinvar_feature_table_dir=$ml_data_dir/clinvar_feature_tables


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
		printf "\n>>> $class:"
		genomic_class_out_file="${gwrvis_distr_dir}/gwrvis.${class}.mutually_excl.bed"

		printf "\n- Getting full genomic class BED file across all chromosomes ... (mutually exclusive with other classes):\n$genomic_class_out_file"
		cat ${gwrvis_distr_dir}/BED/gwrvis_scores_chr*.genomic_coords.${class}.bed | sed "s/\$/\t$class/g" > $genomic_class_out_file
		printf "\n[$class] Total nt length: "`cat $genomic_class_out_file | awk '{sum+=$3-$2} END {print sum}'`

		# Keep track of interval ranges per genomic class
		cat $genomic_class_out_file | awk '{print $3-$2"\t"$1"\t"$2"\t"$3}' | sort -V | awk '{print $2"\t"$3"\t"$4"\t"$1}' > ${gwrvis_distr_dir}/interval_ranges.${class}.mutually_excl.bed


		printf "\n- Retrieving feature table for $class ..."
		full_feature_table_by_genomic_class=${out_feature_table_dir}/full_gwrvis_and_regulatory_features.${class}.tsv

		printf "\n- Keeping coordinates from the genomic class regions (not the full 3kb windows) -- also removing regions with 'NaN' gwRVIS value"
		tail -n +2 $full_feature_table | intersectBed -wao -a $genomic_class_out_file -b stdin| grep $class | grep -v "NaN" | sortBed | cut --complement -f6,7,8,9,37 > ${full_feature_table_by_genomic_class}.tmp


		printf "\nAdding header to feature table by class and saving into file"
		add_header_to_bed_by_genomic_class $full_feature_table_by_genomic_class
		printf "\nOutput file: $full_feature_table_by_genomic_class\n\n"
		
	done
}


function merge_feature_tables_from_all_classes {

	printf "\nMerging feature tables from all genomic classes ..."
	rm -f $out_full_feature_table 

	# Minor bug: remove a single overlapping small region (between an intron and a lincrna) -- calling awk !seen for that ...
	cat ${out_feature_table_dir}/full_gwrvis_and_regulatory_features.*.tsv | grep -v gwrvis | awk '!seen[$1"_"$2] {print} {++seen[$1"_"$2]}' > ${out_full_feature_table}.tmp

	# add header
	cat ${out_feature_table_dir}/full_gwrvis_and_regulatory_features.*.tsv | head -1 > ${out_full_feature_table}.header
	cat ${out_full_feature_table}.header ${out_full_feature_table}.tmp > $out_full_feature_table
	rm ${out_full_feature_table}.tmp ${out_full_feature_table}.header
	printf "\n$out_full_feature_table\n"
}



function add_clinvar_annotation {

	mkdir -p $clinvar_feature_table_dir


	# First, subtract any benign variants from the current pathogenic file!
	subtractBed -a ../other_datasets/variant_annotation/${pathogenic_set}/${pathogenic_set}.pathogenic.bed -b ../other_datasets/variant_annotation/${benign_set}/${benign_set}.benign.bed > $clinvar_feature_table_dir/clean.${pathogenic_set}.pathogenic.bed

	
	# > Get intersections with clinvar pathogenic/benign
	tail -n +2 $out_full_feature_table | intersectBed -wo -a $clinvar_feature_table_dir/clean.${pathogenic_set}.pathogenic.bed -b stdin | cut --complement -f5,6,7,37 > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_pathogenic.bed.tmp
	tail -n +2 $out_full_feature_table | intersectBed -wo -a ../other_datasets/variant_annotation/${benign_set}/${benign_set}.benign.bed -b stdin | cut --complement -f5,6,7,37 > $clinvar_feature_table_dir/full_feature_table.${benign_set}_benign.bed.tmp

	cat $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_pathogenic.bed.tmp $clinvar_feature_table_dir/full_feature_table.${benign_set}_benign.bed.tmp > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed.tmp

	# add header
	tmp_header_file=$clinvar_feature_table_dir/header.tmp
	cat $out_full_feature_table | head -1 | sed 's/gwrvis/clinvar_annot\tgwrvis/' | sed 's/\t$//' > $tmp_header_file
	cat $tmp_header_file $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed.tmp > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed

	# clenaup tmp files
	rm $clinvar_feature_table_dir/*.tmp

	printf "\nOutput file with all annotations: $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed\n"
}



function add_external_genome_wide_scores {

	clinvar_feature_table_bed="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"

	for score in phastCons46way phyloP46way cadd orion ncRVIS; do
		echo $score

		# compile on the fly table with pathogenic and benign variants per score based on the defined pathogenic/benign sets
		cat ../other_datasets/genome-wide-scores/${score}/${score}.${pathogenic_set}_pathogenic.bed ../other_datasets/genome-wide-scores/${score}/${score}.${benign_set}_benign.bed > ../other_datasets/genome-wide-scores/${score}/${score}.${pathogenic_set}_${benign_set}.bed

		tail -n+2 $clinvar_feature_table_bed | intersectBed -wo -a stdin -b ../other_datasets/genome-wide-scores/${score}/${score}.${pathogenic_set}_${benign_set}.bed | cut --complement -f34,35,36,38 > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp
		
		# add header with extra column for the new score
		cat $clinvar_feature_table_bed | head -1 | sed "s/\$/\t$score/" > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header
		cat $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed
		rm $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp
	done

	printf "\nOutput file: $clinvar_feature_table_bed\n"
}



# =============== MAIN RUN ================
printf "\n\n------------\nMerging BED files by genomic class and retrieving respective feature table...\n"
get_feature_table_by_genomic_class

printf "\n\n------------\nMerging feature tables from all genomic classes...\n"
merge_feature_tables_from_all_classes


printf "\n\n------------\nAnnotating full feature table (already with genomic class annotation) with ClinVar pathogenic/bengign variants...\n"
add_clinvar_annotation


printf "\n\n------------\nAdding external genome-wide scores (phastCons46way, phyloP46way, CADD, Orion)...\n"
add_external_genome_wide_scores 
