#!/bin/sh
##SBATCH -J annot
#SBATCH --cpus-per-task=5
#SBATCH --mem=4G
#SBATCH -t 4:00:00


# ==== Read / Infer input arguments ====
config_file=$1
input_classes_file=$2
readarray -t input_classes < $input_classes_file

out_dir=`python custom_utils.py $config_file`
printf "Out dir: ${out_dir}\n"


win_len=`cat $config_file | grep "win_len:" | sed 's/^[ \t]*win_len: //' | sed 's/[ ]*$//'`
pathogenic_set=`cat $config_file | grep "pathogenic_set:" | sed 's/^[ \t]*pathogenic_set: //' | sed 's/[ ]*$//'`
benign_set=`cat $config_file | grep "benign_set:" | sed 's/^[ \t]*benign_set: //' | sed 's/[ ]*$//'`
hg_version=`cat $config_file | grep "hg_version:" | sed 's/^[ \t]*hg_version: //' | sed 's/[ ]*$//'`
labelset_size=`cat $config_file | grep "labelset_size:" | sed 's/^[ \t]*labelset_size: //' | sed 's/[ ]*$//'`
discard_zero_values=`cat $config_file | grep "discard_zero_values:" | sed 's/^[ \t]*discard_zero_values: //' | sed 's/[ ]*$//'`
echo "win_len: *$win_len*"
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


		printf "\n- Retrieving feature table for $class ...\n"
		full_feature_table_by_genomic_class=${out_feature_table_dir}/full_gwrvis_and_regulatory_features.${class}.tsv

		
		printf "\n- Keeping coordinates from the genomic class regions (not the full 3kb windows) -- also removing regions with 'NaN' gwRVIS value"
		tail -n +2 $full_feature_table | intersectBed -wao -a $genomic_class_out_file -b stdin| grep $class | grep -v "NaN" | sortBed | cut --complement -f6,7,8,9,37 > $full_feature_table_by_genomic_class  


		# Input here is: full_gwrvis_and_regulatory_features.ucne.tsv 

		printf "\n- Replacing Ensembl window-based annotation with current BED ranges ..."
		./gwrvis_core/intersect_regulatory_features.sh $full_feature_table_by_genomic_class 


		#printf "\nAdding header to feature table by class and saving into file"
		cp $full_feature_table_by_genomic_class ${full_feature_table_by_genomic_class}.tmp
		add_header_to_bed_by_genomic_class $full_feature_table_by_genomic_class 0 0

		printf "\nOutput file: $full_feature_table_by_genomic_class\n\n"

	done
}


function merge_feature_tables_from_all_classes {

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

	out_file="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"
	printf "\nOutput file with all annotations: $out_file\n"


	# BETA  (enable/disable depending on analysis type: disable when training full model, enable when predicting on a test/validation set to avoid circularity
	# De-duplication: Keep only one entry for windows with multiple variants
	#cp $out_file ${out_file}.with_duplicates_in_windows
	#python deduplicate_overlapping_non_coding_variants.py $out_file 
	#printf "\nOutput file with all annotations (without duplicates): $out_file\n"

}



function add_external_genome_wide_scores {

	clinvar_feature_table_bed="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"

	for score in deepsea eigenPC ncER_10bp cdts linsight phastCons46way phyloP46way cadd dann orion ncRVIS; do
		echo $score

		# compile on the fly table with pathogenic and benign variants per score based on the defined pathogenic/benign sets
		if [ $score != deepsea ]; then  # pre-compiled for DeepSEA
			cat ../other_datasets/genome-wide-scores/${score}/${score}.${pathogenic_set}_pathogenic.bed ../other_datasets/genome-wide-scores/${score}/${score}.${benign_set}_benign.bed > ../other_datasets/genome-wide-scores/${score}/${score}.${pathogenic_set}_${benign_set}.bed
		fi

		tail -n+2 $clinvar_feature_table_bed | intersectBed -wo -a stdin -b ../other_datasets/genome-wide-scores/${score}/${score}.${pathogenic_set}_${benign_set}.bed | cut --complement -f34,35,36,38 > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp
		
		# add header with extra column for the new score
		cat $clinvar_feature_table_bed | head -1 | sed "s/\$/\t$score/" > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header
		cat $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp > $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed
		rm $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp
	done

	#echo $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.header 
	#echo $clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.${score}.bed.tmp 

	printf "\nOutput file: $clinvar_feature_table_bed\n"
}









function update_single_nt_gwrvis_per_chr {

	chr=$1
	clinvar_full_feature_table=$2

	single_nt_base_dir="../out/topmed-single-nt-gwrvis-winlen_${win_len}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/single_nt_gwrvis_per_chr/"
	single_nt_gwrvis_chr_bed="$single_nt_base_dir/gwrvis_single_nt.chr${chr}.bed.gz"

	
	tail -n+2 $clinvar_full_feature_table | awk -F '\t' -v ch="$chr" '{if($1 == "chr"ch) print $1"\t"$2"\t"$3}' > "${clinvar_full_feature_table}.chr${chr}.bed"

	tabix $single_nt_gwrvis_chr_bed -B "${clinvar_full_feature_table}.chr${chr}.bed" > "${clinvar_full_feature_table}.chr${chr}.single_nt_gwrvis.bed"

	
	# cleanup
	rm -f ${clinvar_full_feature_table}.chr${chr}.bed

	echo "Done (single_nt_gwrvis) - chr: $chr"

}




function update_single_nt_gwrvis {

	clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"


	for chr in `seq 1 22`; do
	#for chr in 21; do
		echo "Chr: $chr"
		update_single_nt_gwrvis_per_chr $chr $clinvar_full_feature_table &
	done
	wait


	cat  ${clinvar_full_feature_table}.chr*.single_nt_gwrvis.bed > "${clinvar_full_feature_table}.all_chr.single_nt_gwrvis.bed"
	echo "${clinvar_full_feature_table}.all_chr.single_nt_gwrvis.bed"

	
}







function subset_phastcons_primate {

	chr=$1
	clinvar_full_feature_table=$2
	phastCons_primates_chr_bed="../other_datasets/conservation/phastCons46way_primates/bed/chr${chr}.phastCons46way.primates.high_conf_regions.bed.gz"

	tail -n+2 $clinvar_full_feature_table | awk -F '\t' -v ch="$chr" '{if($1 == "chr"ch) print $1"\t"$2"\t"$3}' > "${clinvar_full_feature_table}.chr${chr}.bed"

	tabix $phastCons_primates_chr_bed -B "${clinvar_full_feature_table}.chr${chr}.bed" > "${clinvar_full_feature_table}.chr${chr}.phastcons_primate.bed"

	echo "${clinvar_full_feature_table}.chr${chr}.phastcons_primate.bed"
	

	# TODO: Keep NAs to later impute with median
	tail -n+2 $clinvar_full_feature_table | awk -v ch="$chr" '{if($1 == "chr"ch) print $0}' | intersectBed -wao -a stdin -b "${clinvar_full_feature_table}.chr${chr}.phastcons_primate.bed" | cut -f1-33,37 > "${clinvar_full_feature_table}.chr${chr}.bed"


	# cleanup
	rm "${clinvar_full_feature_table}.chr${chr}.phastcons_primate.bed" 
	echo "[DONE] subset_phastcons_primate - chr $chr"
}




function add_conservation_features {

	clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"

	cat $clinvar_full_feature_table | head -1 > ${clinvar_full_feature_table}.header
	sed -i 's/$/\tphastCons_primate/' ${clinvar_full_feature_table}.header


	for chr in `seq 1 22`; do
	#for chr in 22; do
		echo "Chr: $chr"
		subset_phastcons_primate $chr $clinvar_full_feature_table &
	done
	wait


	cat ${clinvar_full_feature_table}.chr*.bed > "${clinvar_full_feature_table}.all_chr.bed"
	#echo "${clinvar_full_feature_table}.all_chr.bed"

	cat ${clinvar_full_feature_table}.header "${clinvar_full_feature_table}.all_chr.bed" > "${clinvar_full_feature_table}.full_with_phastcons.bed"
	echo -e "Complete.\n${clinvar_full_feature_table}.full_with_phastcons.bed"
	

	#cp $clinvar_full_feature_table ${clinvar_full_feature_table}.original
	cp ${clinvar_full_feature_table}.full_with_phastcons.bed $clinvar_full_feature_table
		
}



function add_distance_from_closest_tss {
	
	clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"
	tail -n+2 $clinvar_full_feature_table | cut -f1,2,3 > ${clinvar_full_feature_table}.tmp

	tss_ref_file="../ensembl/ensembl_genes.hg19.TSS.bed"


	echo "TSS_distance" > ${clinvar_full_feature_table}.tmp.tss
	closest-features --closest --delim '\t' ${clinvar_full_feature_table}.tmp $tss_ref_file | awk 'function abs(v) {return v < 0 ? -v : v} {print abs($5-$2)}' >> ${clinvar_full_feature_table}.tmp.tss

	

	paste $clinvar_full_feature_table ${clinvar_full_feature_table}.tmp.tss > ${clinvar_full_feature_table}.tmp

	mv ${clinvar_full_feature_table}.tmp $clinvar_full_feature_table
	echo $clinvar_full_feature_table


	#cleanup
	rm ${clinvar_full_feature_table}.tmp.tss

}











# =============== MAIN RUN ================
printf "\n\n------------\nMerging BED files by genomic class and retrieving respective feature table...\n"
get_feature_table_by_genomic_class



printf "\n\n------------\nMerging feature tables from all genomic classes...\n"
merge_feature_tables_from_all_classes


printf "\n\n------------\nAnnotating full feature table (already with genomic class annotation) with ClinVar pathogenic/bengign variants...\n"
add_clinvar_annotation




clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"
# -- REDUNDANT
#printf "\n\n------------\nGetting single-nt Ensembl annotation..."
#./gwrvis_core/full_table_intersect_regulatory_features.sh $clinvar_full_feature_table



printf "\n\n------------\nAdding external genome-wide scores (phastCons46way, phyloP46way, CADD, Orion)...\n"
add_external_genome_wide_scores 




printf "\n\n------------\nAdding primate-phastCons annotation...\n"
add_conservation_features


printf "\n\n-----------\nAdd distance of each variant from closest TSS as feature...\n"
add_distance_from_closest_tss




printf "\n\n-----------\nSort clinvar full feature table (patho_benign) set by chr and start coord...\n"
cat $clinvar_full_feature_table | head -1 > ${clinvar_full_feature_table}.tmp
tail -n+2 $clinvar_full_feature_table | sort -k1,1V -k2,2n >> ${clinvar_full_feature_table}.tmp
mv ${clinvar_full_feature_table}.tmp $clinvar_full_feature_table



# -- REDUNDANT
#printf "\n\n-----------\nIntersect variants with genomic classes on single-nt level (by-passing window-based gwRVIS approach)...\n"
#gwrvis_core/intersect_variant_table_with_genomic_classes.sh $config_file $input_classes_file



:'
printf "\n\n-----------\nGet single-nt gwRVIS values from all intervals in the clinvar_full_feature table...\n"
update_single_nt_gwrvis



printf "\n\n-----------\nUpdate clinvar_feature_table with single-nt resolution gwrvis scores...\n"
python -u gwrvis_core/update_clinvar_table_with_singlent_gwrvis.py $config_file
'
