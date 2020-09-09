#!/bin/sh

in_file="ensembl_genes.hg19.gz"
out_full_file="ensembl_genes.hg19.TSS_genenames.txt"
out_file="ensembl_genes.hg19.TSS.bed"


zcat < $in_file | tail -n+2 | awk ' BEGIN{OFS="\t"}{ if($4=="+") {print $3,$5,$5+1,$2 "_" $13,".",$4}  else {print $3,$6-1,$6,$2 "_" $13,".",$4}  }' > $out_full_file


cat $out_full_file | cut -f1,2,3 | sortBed | mergeBed > $out_file



