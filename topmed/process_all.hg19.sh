#!/bin/bash -l
#SBATCH -J process_topmed_hg19  # Set name for current job 
#SBATCH -o out.process_topmed_hg19  # Set job output log 
#SBATCH -e err.process_topmed_hg19  # Set job error log 
#SBATCH --cpus-per-task=10         # Request 9 CPUs (cores) on a single node 
#SBATCH --mem=24000          # Request amount of memory 
#SBATCH -t 24:0:0            # Request 24 hours runtime

FILTER=$1
SNV_ONLY=$2

if [ $# -ne 2 ]; then         
	echo "[Error]: incorrect number of input arguments."
	echo ">> Expected call format: ./process_all.hg19.sh [FILTER: 0|1] [SNV_ONLY: 0|1]"         
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

	echo "(1/5) Retrieving CHR, POS, REF, ALT, FILTER ..."
	cat vcf/bravo-dbsnp-all.hg19.vcf | grep -v '##' | grep -E 'PASS|^#' | cut -f1,2,4,5,7 > topmed_subtable${filtered_annot}.txt
	echo "(2/5) Retrieving INFO fields (AC, AN, AF, DP)..."
	# Appending a placeholder value=500000 for the DP column in the end, for compatibility with processing of VCF files for gwRVIS calculation
	cat vcf/bravo-dbsnp-all.hg19.vcf | grep -v '##' | grep -E 'PASS|^#' | cut -f8 | awk -F ";|=" '{print $2"\t"$4"\t"$6"\t"500000}' > topmed_info${filtered_annot}.txt
else
	echo ">   [No filter applied; retaining all entries]"

	echo "(1/5) Retrieving CHR, POS, REF, ALT, FILTER ..."
	cat vcf/bravo-dbsnp-all.hg19.vcf | grep -v '##' | cut -f1,2,4,5,7 > topmed_subtable${filtered_annot}.txt
	echo "(2/5) Retrieving INFO fields (AC, AN, AF, DP)..."
	cat vcf/bravo-dbsnp-all.hg19.vcf | grep -v '##' | cut -f8 | awk -F ";|=" '{print $2"\t"$4"\t"$6"\t"500000}' > topmed_info${filtered_annot}.txt

fi


out_table_name=topmed_table${filtered_annot}
out_table_vcf=${out_table_name}.vcf
echo "(3/5) Combinining the two tables ..."
paste topmed_subtable${filtered_annot}.txt topmed_info${filtered_annot}.txt > ${out_table_vcf}.tmp


# Keep only SNVs if defined in input args
if [ "$SNV_ONLY" -eq 1 ]; then
	echo ">   [retaining SNVs only]"
	cat  ${out_table_vcf}.tmp | awk '{ if(NR==1) print $0; else if ( (length($3)==length($4)) && (length($3)==1) ) print $0 }' > ${out_table_vcf}
	rm ${out_table_vcf}.tmp
else
	mv ${out_table_vcf}.tmp ${out_table_vcf}
fi


# [Deprecated] Convert topmed table from VCF to BED: requires an extra column with coord vcf_coord-1, 
# to go from 1-based to 0-based coordinate system
#out_table_bed=${out_table_name}.bed
#cat $out_table_vcf | awk '{if($1 ~ /^#/) print $0; else print $1"\t"$2-1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' > $out_table_bed

## Liftover using Picard from gatk: from hg38 to hg37/19
#echo "(4/5) Liftover from hg38 to hg19 - using gatk ..."
#./gatk_liftover.sh


echo "(5/5) Splitting VCF entries by chromosome ..."
python split_table_into_chr.py $out_table_vcf 1 1	# args: topmed_processed_table_file, FILTER [0/1], SNVs_ONLY [0/1]


#echo "(5/5) Cleaning up temporary files ..."
rm topmed_subtable${filtered_annot}.txt topmed_info${filtered_annot}.txt 
#rm topmed_table.txt
