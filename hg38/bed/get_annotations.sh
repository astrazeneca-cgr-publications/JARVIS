#!/bin/bash

# Author: Dimitrios Vitsios
# Based on: https://davetang.org/muse/2013/01/18/defining-genomic-regions/


v=29  # ensembl release version
hg_v=38	# human genome assembly version


# Download GTF annotation file
gtf_file=gencode.v${v}.annotation.gtf.gz
if [ ! -f $gtf_file ]; then
	echo "Downloading gtf file: gencode.v${v}.annotation.gtf.gz"
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${v}/gencode.v${v}.annotation.gtf.gz
fi


#echo "GTF summary:"
#zcat gencode.v29.annotation.gtf.gz | grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn

cds_file=gencode_v${v}_cds_merged.bed
utr_file=gencode_v${v}_utr_merged.bed
exon_file=gencode_v${v}_exon_merged.bed
intron_file=gencode_v${v}_intron.bed
intergenic_file=gencode_v${v}_intergenic.bed


# ----  CCDS  ----
if [ ! -f $cds_file ]; then
	echo "Extracting CDS file ..."
	zcat gencode.v${v}.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="CDS" {print $1,$4-1,$5}' | sortBed | mergeBed -i - > $cds_file
fi


# ----  UTR  ----
if [ ! -f $utr_file ]; then
	echo "Extracting UTR file ..."
	zcat gencode.v${v}.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="UTR" {print $1,$4-1,$5}' | sortBed | mergeBed -i - > $utr_file
fi


# ----  Exons  ---- (intermediate file)
if [ ! -f $exon_file ]; then
	echo "Extracting exons file ..."
	zcat gencode.v${v}.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sortBed | mergeBed -i - > $exon_file
fi

# ---- Introns  ----
if [ ! -f $intron_file ]; then
	echo "Extracting introns file ..."
	zcat gencode.v${v}.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | subtractBed -a stdin -b $exon_file > $intron_file
fi

echo "Check intersection between exon and intron file (should be empty):"
intersectBed -a $exon_file -b $intron_file
echo "... OK"


# ----  Intergenic regions  ----
if [ ! -f hg${hg_v}.genome ]; then
	echo "Extracting hg${hg_v}.genome file ..."
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg${hg_v}.chromInfo"  > hg${hg_v}.genome
fi
#echo "\n[Important]: fix manually the order of chromosomes in hg${hg_v}.genome file, based on the chrom order in gtf after sortBed"
#echo -e "The correct order can be extracted by running:\nzcat gencode.v${v}.annotation.gtf.gz | awk 'BEGIN{OFS=\"\t\";} $3==\"gene\" {print $1,$4-1,$5}' | sortBed | cut -f1 | uniq"

if [ ! -f $intergenic_file ]; then
	echo "Extracting intergenic file ..."

	zcat gencode.v${v}.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | complementBed -i stdin -g hg${hg_v}.genome > $intergenic_file

fi


