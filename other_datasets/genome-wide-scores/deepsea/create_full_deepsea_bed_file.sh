#!/bin/bash

base_dir="genomic-files/deepsea-vcf"


for dataset in clinvar_denovodb ncER-generalization_other_denovodb ncER-generalization_ncRNA_denovodb ncER-mendelian_denovodb ncER-gwas_catalog_denovodb; do

	cur_dir="${base_dir}/${dataset}-parts"

	out_file="deepsea.${dataset}.bed"
	rm -f $out_file

	for f in `ls ${cur_dir}/infile.vcf.out.snpclass*`; do
	
		echo $f

		tail -n+2 $f | awk -F"," '{print $2"\t"$3-1"\t"$3"\t"$9}' >> $out_file
	
	done

	echo "Out file: $out_file"

done

