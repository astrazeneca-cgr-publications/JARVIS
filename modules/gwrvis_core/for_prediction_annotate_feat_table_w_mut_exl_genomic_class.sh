#!/bin/sh
#SBATCH -J for_pred_annot
#SBATCH -o for_pred_annot.out
#SBATCH --cpus-per-task=5
#SBATCH --mem=200G
#SBATCH -t 24:00:00


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

		# Replace Ensembl window-based annotation with current BED ranges
		./gwrvis_core/intersect_regulatory_features.sh $full_feature_table_by_genomic_class


		printf "\n- Keeping coordinates from the genomic class regions (not the full 3kb windows) -- also removing regions with 'NaN' gwRVIS value"
		#tail -n +2 $full_feature_table | intersectBed -wao -a $genomic_class_out_file -b stdin| grep $class | grep -v "NaN" | sortBed | cut --complement -f6,7,8,9,37 > ${full_feature_table_by_genomic_class}.tmp
		#echo ${full_feature_table_by_genomic_class}.tmp


		#printf "\nAdding header to feature table by class and saving into file"
		#add_header_to_bed_by_genomic_class $full_feature_table_by_genomic_class 0 0
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








function add_mock_clinvar_annotation {



	out_file="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"

	python gwrvis_core/for_testing_prepare_genome_wide_mock_clinvar_table.py $out_full_feature_table $out_file
	printf "\nOutput file with all annotations: $out_file\n"

	#python gwrvis_core/for_testing_prepare_single_nt_genome_wide_bed.py $out_full_feature_table $out_feature_table_dir	
}










function update_single_nt_gwrvis_per_chr {

	chr=$1
	clinvar_full_feature_table=$2

	single_nt_base_dir="../out/topmed-single-nt-gwrvis-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/single_nt_gwrvis_per_chr/"
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

}









function update_single_nt_ensembl_annot_per_chr {

	chr=$1
	clinvar_full_feature_table=$2


	single_nt_chr_bed="${clinvar_full_feature_table}.chr${chr}.single_nt_gwrvis.bed"
	# -- DEBUG: single_nt_chr_bed="chr21.tmp"
	

	# reset to original bed:
	cat $single_nt_chr_bed | cut -f1-4 > ${single_nt_chr_bed}.tmp
	mv ${single_nt_chr_bed}.tmp $single_nt_chr_bed




	cell_line="Monocytes_CD14plus"
	for regul_elem in CTCF_Binding_Site Enhancer Open_chromatin TF_binding_site H3K27ac H3K27me3 H4K20me1 H3K9ac H3K4me1 H3K4me2 H3K4me3 H3K36me3; do   	

		echo $regul_elem

		ensembl_regul_file="../ensembl/GRCh37-Regulatory_Features/${cell_line}.${regul_elem}.sorted.bed"
		echo "ensembl_regul_file: $ensembl_regul_file"

		tmp_regul_inters_bed=${single_nt_chr_bed}_${regul_elem}.tmp



		
		intersectBed -wao -a $single_nt_chr_bed -b $ensembl_regul_file |  awk '!seen[$1"_"$2] {print} {++seen[$1"_"$2]}' | rev | cut -f1 | rev | awk '{if($1 > 0) {print 1} else print 0}' > $tmp_regul_inters_bed

		echo "tmp_regul_inters_bed:"
		echo `cat $tmp_regul_inters_bed | wc -l `


		paste $single_nt_chr_bed $tmp_regul_inters_bed > ${single_nt_chr_bed}.tmp
		mv ${single_nt_chr_bed}.tmp $single_nt_chr_bed


		echo -e "\nDone - chr$chr, $regul_elem"
		echo `cat $single_nt_chr_bed | wc -l`

		rm $tmp_regul_inters_bed



	done



	echo -e "chr\tstart\tend\tgwrvis\tCTCF_Binding_Site\tEnhancer\tOpen_chromatin\tTF_binding_site\tH3K27ac\tH3K27me3\tH4K20me1\tH3K9ac\tH3K4me1\tH3K4me2\tH3K4me3\tH3K36me3" > ${single_nt_chr_bed}.header	

	# add header
	cat ${single_nt_chr_bed}.header $single_nt_chr_bed > ${single_nt_chr_bed}.tmp
	mv ${single_nt_chr_bed}.tmp $single_nt_chr_bed


	# cleanup
	rm ${single_nt_chr_bed}.header

}



function update_single_nt_ensembl_annot {

	clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"


	cnt=0
	for chr in `seq 1 22`; do 
	#for chr in `seq 21 21`; do 
		echo "Chr: $chr"
		update_single_nt_ensembl_annot_per_chr $chr $clinvar_full_feature_table &
	


		if [ "$cnt" == 5 ]; then
			cnt=0
			wait
		fi

		cnt=$((cnt+1))

	done
	wait

}






# =============== MAIN RUN ================
printf "\n\n------------\nMerging BED files by genomic class and retrieving respective feature table...\n"
#get_feature_table_by_genomic_class



printf "\n\n------------\nMerging feature tables from all genomic classes...\n"
#merge_feature_tables_from_all_classes


printf "\n\n------------\nAnnotating full feature table (already with genomic class annotation) with ClinVar pathogenic/bengign variants...\n"
# add_mock_clinvar_annotation





# {Add primate-phastCons annotation
# Add distance of each variant from closest TSS as feature}






printf "\n\n------------\nGet single-nt gwRVIS values from all intervals in the clinvar_full_feature table...\n"
update_single_nt_gwrvis





printf "\n\n------------\nGet single-nt Ensembl annotation...\n"
update_single_nt_ensembl_annot







# NEW - this is a critical step for getting single-nt Ensembl annotation
# {{clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed" 
# {{./gwrvis_core/full_table_intersect_regulatory_features.sh $clinvar_full_feature_table







