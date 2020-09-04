#!/bin/sh

# ================ Read / Infer input arguments ================
config_file=$1
input_classes_file=$2
readarray -t input_classes < $input_classes_file

out_dir=`python custom_utils.py $config_file`
printf "Out dir: ${out_dir}\n"


pathogenic_set=`cat $config_file | grep "pathogenic_set:" | sed 's/^[ \t]*pathogenic_set: //' | sed 's/[ ]*$//'`
benign_set=`cat $config_file | grep "benign_set:" | sed 's/^[ \t]*benign_set: //' | sed 's/[ ]*$//'`
echo "Pathogenic set: *$pathogenic_set*"
echo "Benign set: *$benign_set*"


# ==== Global variables ====
ml_data_dir="${out_dir}/ml_data"
mkdir -p  $ml_data_dir

out_feature_table_dir=$ml_data_dir/feature_tables
mkdir -p $out_feature_table_dir

clinvar_feature_table_dir=$ml_data_dir/clinvar_feature_tables
mkdir -p $clinvar_feature_table_dir









clinvar_full_feature_table="$clinvar_feature_table_dir/full_feature_table.${pathogenic_set}_${benign_set}.bed"



remaining_clinvar_table="${clinvar_full_feature_table}.tmp"
tail -n+2 $clinvar_full_feature_table > $remaining_clinvar_table


declare -A bed_files
while IFS=$'\t' read -r key value priority; do
	bed_files[$key]=$value
done < static_files/genomic_classes.log


echo "${bed_files['ccds']}"


for class in "${input_classes[@]}"; do

	intersectBed -a $remaining_clinvar_table -b ${bed_files[$class]} > ${remaining_clinvar_table}.${class}_tmp
	subtractBed -a $remaining_clinvar_table -b ${remaining_clinvar_table}.${class}_tmp > ${remaining_clinvar_table}.tmp

	mv ${remaining_clinvar_table}.tmp $remaining_clinvar_table

done

# cleanup 
rm ${remaining_clinvar_table}.tmp


echo ${remaining_clinvar_table}.${class}_tmp



# Call python module to update genomic classes in the clinvar_full_feature_table
python gwrvis_core/get_genomic_classes_for_jarvis.py $config_file $input_classes_file
