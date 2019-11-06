#!/bin/bash -l
#SBATCH -J process_topmed_hg38  # Set name for current job 
#SBATCH -o out.process_topmed_hg38  # Set job output log 
#SBATCH -e err.process_topmed_hg38  # Set job error log 
#SBATCH --cpus-per-task=1         # Request 9 CPUs (cores) on a single node 
#SBATCH --mem=24G          # Request amount of memory 
#SBATCH -t 24:0:0            # Request 24 hours runtime

FILTER=$1
SNV_ONLY=$2

if [ $# -ne 2 ]; then         
	echo "[Error]: incorrect number of input arguments."
	echo ">> Expected call format: ./process_all.sh [FILTER: 0|1] [SNV_ONLY: 0|1]"         
	exit 
fi


filtered_annot=''
if [ "$FILTER" -eq 1 ]; then
	filtered_annot='.filtered'
fi
if [ "$SNV_ONLY" -eq 1 ]; then
	filtered_annot=".SNVs${filtered_annot}"
fi


if [ "$FILTER" -eq 1 ]; then
	echo ">   [retaining only entries with PASS flag]"

	echo "(1/4) Retrieving CHR, POS, REF, ALT, FILTER ..."
	cat vcf/bravo-dbsnp-all.vcf | grep -v '##' | grep -E 'PASS|^#' | cut -f1,2,4,5,7 > topmed_subtable${filtered_annot}.txt
	echo "(2/4) Retrieving INFO fields (AC, AN, AF, DP)..."
	# Appending a placeholder value=500000 for the DP column in the end, for compatibility with processing of VCF files for gwRVIS calculation
	cat vcf/bravo-dbsnp-all.vcf | grep -v '##' | grep -E 'PASS|^#' | cut -f8 | awk -F ";|=" '{print $8"\t"$10"\t"$6"\t"500000}' > topmed_info${filtered_annot}.txt
else
	echo ">   [No filter applied; retaining all entries]"

	echo "(1/4) Retrieving CHR, POS, REF, ALT, FILTER ..."
	cat vcf/bravo-dbsnp-all.vcf | grep -v '##' | cut -f1,2,4,5,7 > topmed_subtable${filtered_annot}.txt
	echo "(2/4) Retrieving INFO fields (AC, AN, AF, DP)..."
	cat vcf/bravo-dbsnp-all.vcf | grep -v '##' | cut -f8 | awk -F ";|=" '{print $8"\t"$10"\t"$6"\t"500000}' > topmed_info${filtered_annot}.txt

fi


out_table_name=topmed_table${filtered_annot}
out_table_vcf=${out_table_name}.vcf
echo "(3/4) Combinining the two tables ..."
paste topmed_subtable${filtered_annot}.txt topmed_info${filtered_annot}.txt > ${out_table_vcf}.tmp


# Keep only SNVs if defined in input args
if [ "$SNV_ONLY" -eq 1 ]; then
	echo ">   [retaining SNVs only]"
	cat  ${out_table_vcf}.tmp | awk '{ if(NR==1) print $0; else if ( (length($3)==length($4)) && (length($3)==1) ) print $0 }' > ${out_table_vcf}
	rm ${out_table_vcf}.tmp
else
	mv ${out_table_vcf}.tmp ${out_table_vcf}
fi



echo "(4/4) Splitting VCF entries by chromosome ..."
python split_table_into_chr.py $out_table_vcf $FILTER $SNV_ONLY	# args: topmed_processed_table_file, FILTER [0/1], SNVs_ONLY [0/1]


#echo "(5/5) Cleaning up temporary files ..."
rm topmed_subtable${filtered_annot}.txt topmed_info${filtered_annot}.txt 
#rm topmed_table.txt
