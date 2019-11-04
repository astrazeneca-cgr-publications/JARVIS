#!/bin/bash

cell_line="Monocytes_CD14plus"

echo "TF"
tf_file="Monocytes_CD14plus.TF_binding_sites.martquery_0702204213_524.txt.gz"
zcat < $tf_file | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\tTF_binding_site"}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.TF_binding_sites.sorted.bed



echo "Open chromatin"
open_chromatin_file="Monocytes_CD14plus.Open_chromatin.txt.gz"
zcat < $open_chromatin_file | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\tOpen_chromatin"}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.Open_chromatin.sorted.bed



echo "Enhancer"
enhancer_file="Monocytes_CD14plus.Enhancers.martquery_0702203843_872.txt.gz"
zcat < $enhancer_file | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\tEnhancer"}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.Enhancers.sorted.bed



echo "CTCF"
ctcf_file="Monocytes_CD14plus.CTCF_binding_sites.martquery_0702204031_788.txt.gz"
zcat < $ctcf_file | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\tCTCF_Binding_Site"}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.CTCF_binding_sites.sorted.bed



echo "H3K27ac"
H3K27ac_file="Monocytes_CD14plus.H3K27ac.mart_export.txt.gz"
zcat < $H3K27ac_file | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\tH3K27ac"}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.H3K27ac.sorted.bed



echo "H3K27me3"
H3K27me3_file="Monocytes_CD14plus.H3K27me3.mart_export.txt.gz"
zcat < $H3K27me3_file | tail -n+2 | awk '{print $1"\t"$2"\t"$3"\tH3K27me3"}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.H3K27me3.sorted.bed



echo "Histone Modifications:"
histone_mod_file="Monocytes_CD14plus.Histone_modifications.mart_export.txt.gz"
for histone_mod in H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K9ac H4K20me1; do
	echo "- $histone_mod"

	zcat < $histone_mod_file | tail -n+2 | grep $histone_mod | awk -v hist_mod="$histone_mod" '{print $1"\t"$2"\t"$3"\t"hist_mod}' | sortBed | mergeBed -c 4 -o 'distinct' | sed 's/^/chr/g' > ${cell_line}.${histone_mod}.sorted.bed

done

