#!/bin/sh

# ==== Read / Infer input arguments ====
config_file=$1
input_classes_file=$2
readarray -t input_classes < $input_classes_file

out_dir=`python custom_utils.py $config_file`
printf "Out dir: ${out_dir}\n"


pathogenic_set=`cat $config_file | grep "pathogenic_set:" | sed 's/^[ \t]*pathogenic_set: //' | sed 's/[ ]*$//'`
benign_set=`cat $config_file | grep "benign_set:" | sed 's/^[ \t]*benign_set: //' | sed 's/[ ]*$//'`
hg_version=`cat $config_file | grep "hg_version:" | sed 's/^[ \t]*hg_version: //' | sed 's/[ ]*$//'`
labelset_size=`cat $config_file | grep "labelset_size:" | sed 's/^[ \t]*labelset_size: //' | sed 's/[ ]*$//'`
discard_zero_values=`cat $config_file | grep "discard_zero_values:" | sed 's/^[ \t]*discard_zero_values: //' | sed 's/[ ]*$//'`
echo "Pathogenic set: *$pathogenic_set*"
echo "Benign set: *$benign_set*"
echo "hg version: *$hg_version*"
echo "labelset_size: *$labelset_size*"
echo "discard_zero_values: *$discard_zero_values*"


variant_annot_dir="../other_datasets/variant_annotation-${hg_version}"


# ==== Global variables ====
gwrvis_distr_dir="${out_dir}/gwrvis_distribution"
full_feature_table="${out_dir}/full_genome_out/full_gwrvis_and_regulatory_features.bed"

ml_data_dir="${out_dir}/ml_data"
mkdir -p  $ml_data_dir

out_feature_table_dir=$ml_data_dir/feature_tables
mkdir -p $out_feature_table_dir

clinvar_feature_table_dir=$ml_data_dir/clinvar_feature_tables
mkdir -p $clinvar_feature_table_dir
conservation_feature_table_dir=$ml_data_dir/conservation_feature_tables
mkdir -p $conservation_feature_table_dir
tmp_conservation_out=$conservation_feature_table_dir/tmp
mkdir -p $tmp_conservation_out

# Output fle
out_full_feature_table="${out_feature_table_dir}/full_gwrvis_and_regulatory_features.All_genomic_classes.tsv"




function add_header_to_bed_by_genomic_class () {

	full_out_table=$1
	is_conservation_table=$2
	append_score_colname=$3

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

	printf "%s\t" "${full_header[@]}" > ${full_out_table}.header

	if [ "$is_conservation_table" -eq "1" ]; then
		sed -i 's/gwrvis/conservation\_annot\tgwrvis/' ${full_out_table}.header
		sed -i 's/\t$//' ${full_out_table}.header
	fi
	if [ "$append_score_colname" -eq "1" ]; then
		sed -i "s/$/\t${score}/" ${full_out_table}.header
	fi

	sed -i 's/$/\n/' ${full_out_table}.header
	echo -e "\n${full_out_table}.header"

	cat ${full_out_table}.header ${full_out_table}.tmp > ${full_out_table}
	rm ${full_out_table}.header ${full_out_table}.tmp
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
		add_header_to_bed_by_genomic_class $full_feature_table_by_genomic_class 0 0
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





function add_conservation_annotation {


	# Add D{labelset_size}[.no_zeros] annotation -- read it from config
	file_annot="D${labelset_size}" 
	if [ $discard_zero_values == "1" ]; then
		file_annot="${file_annot}.no_zeros"
	fi
	primate_conservation_dir="../other_datasets/conservation/phastCons46way_primates/most_least_conserved_files_by_class/${file_annot}"
		
	for genomic_class in "${input_classes[@]}"; do
	#for genomic_class in intron; do
		echo -e "\n> ${genomic_class}:"

		most_conserved_file="$primate_conservation_dir/${genomic_class}.most_conserved.${file_annot}.bed"
		least_conserved_file="$primate_conservation_dir/${genomic_class}.least_conserved.${file_annot}.bed"

		# Get regions with ambiguous conservation profile (i.e. both conserved and non-conserved regions intersecting the same gwRVIS regions)
		echo "Extracting regions with ambiguous conservation profile ..."
		full_table="${out_feature_table_dir}/full_gwrvis_and_regulatory_features.${genomic_class}.tsv"
		full_table_most_conserved="${tmp_conservation_out}/full_table.${genomic_class}.most_conserved.${file_annot}.bed"
		full_table_least_conserved="${tmp_conservation_out}/full_table.${genomic_class}.least_conserved.${file_annot}.bed"

		tail -n+2 $full_table | intersectBed -wa -a stdin -b $most_conserved_file > $full_table_most_conserved
		tail -n+2 $full_table | intersectBed -wa -a stdin -b $least_conserved_file > $full_table_least_conserved

		ambiguous_conserv_file="${tmp_conservation_out}/${genomic_class}.ambiguous_conservation.${file_annot}.bed"
		intersectBed -a $full_table_most_conserved -b $full_table_least_conserved > $ambiguous_conserv_file


		# Create conservation annotation file for current genomic class -- without ambiguous conservation regions
		echo "Annotate full feature table with unambiguous conservation labels ..."
		concat_label_file="${tmp_conservation_out}/${genomic_class}.conservation_labels.${file_annot}.bed"
		cat $least_conserved_file $most_conserved_file | subtractBed -a stdin -b $ambiguous_conserv_file > $concat_label_file
		#rm $least_conserved_file $most_conserved_file
		echo "Output file: $concat_label_file"


		# Annotate full feature table with unambiguous conservation annotations
		echo "Intersecting with full feature table for $genomic_class class ..."
		tail -n+2 $full_table | intersectBed -wb -a $concat_label_file -b stdin | cut --complement -f5,6,7 > $conservation_feature_table_dir/full_feature_table.conservation.${genomic_class}.${file_annot}.bed 

	done


	full_out_file="$conservation_feature_table_dir/full_feature_table.conservation.All_genomic_classes.${file_annot}.bed"

	# Cleanup previously created files
	rm -f $full_out_file
	for score in orion; do
		score_out_file="$conservation_feature_table_dir/full_feature_table.conservation.All_genomic_classes.${score}.${file_annot}.bed"
		rm -f $score_out_file
	done
	cat $conservation_feature_table_dir/full_feature_table.conservation.*.${file_annot}.bed > ${full_out_file}.tmp


	add_header_to_bed_by_genomic_class $full_out_file 1 0
	printf "\nOutput file: $full_out_file\n\n"


	# Get conservation-annotated tables for other genome-wide scores
	for score in orion; do
		echo -e "\nScore: $score"	
		score_conservation_ref_file="../other_datasets/genome-wide-scores/${score}/${score}.conservation_annotation.${file_annot}.bed"
		score_out_file="$conservation_feature_table_dir/full_feature_table.conservation.All_genomic_classes.${score}.${file_annot}.bed"

		tail -n+2 $full_out_file | intersectBed -wo -a stdin -b $score_conservation_ref_file | rev | cut --complement -f1,3,4,5 | rev > ${score_out_file}.tmp
	
		add_header_to_bed_by_genomic_class $score_out_file 1 1

		echo "Out file: $score_out_file"
	done

}






function add_clinvar_annotation {

	# First, subtract any benign variants from the current pathogenic file!
	subtractBed -a ${variant_annot_dir}/${pathogenic_set}/${pathogenic_set}.pathogenic.bed -b ${variant_annot_dir}/${benign_set}/${benign_set}.benign.bed > $clinvar_feature_table_dir/clean.${pathogenic_set}.pathogenic.bed

	
	# > Get intersections with clinvar pathogenic/benign
	tail -n +2 $out_full_feature_table | intersectBed -wo -a $clinvar_feature_table_dir/clean.${pathogenic_set}.pathogenic.bed -b stdin | cut --complement -f5,6,7,37 > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_pathogenic.bed.tmp
	tail -n +2 $out_full_feature_table | intersectBed -wo -a ${variant_annot_dir}/${benign_set}/${benign_set}.benign.bed -b stdin | cut --complement -f5,6,7,37 > $clinvar_feature_table_dir/full_feature_table.${benign_set}_benign.bed.tmp

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

	for score in phastCons46way phyloP46way cadd dann orion ncRVIS; do
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

#add_conservation_annotation

printf "\n\n------------\nAnnotating full feature table (already with genomic class annotation) with ClinVar pathogenic/bengign variants...\n"
add_clinvar_annotation


printf "\n\n------------\nAdding external genome-wide scores (phastCons46way, phyloP46way, CADD, Orion)...\n"
add_external_genome_wide_scores 
