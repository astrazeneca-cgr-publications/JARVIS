#!/bin/bash

vcf_dir=$1
chr=$2
KEEP_PASS_ONLY=$3
FILTER_SEGDUP=$4
FILTER_LCR=$5
out_dir=$6
population=$7


if [ $# -ne 7 ]; then
	echo "[Error]: incorrect number of input arguments."
	echo ">> Expected call format: ./get_chr_table.sh [vcf_dir] [chr] [KEEP_PASS_ONLY: 0|1] [FILTER_SEGDUP: 0|1] [FILTER_LCR: 0|1] [out_dir] [population]"
	exit
fi
mkdir -p $out_dir
tmp_dir=$out_dir/tmp
mkdir -p $tmp_dir


# >> Define parameters for keeping/filtering certain types of variants based on annotation for PASS, segdup or lcr
filter_params="grep -E ';AF=|^#'"   # exclude entries that don't have AF information (i.e. have AC,AN=0); also keep header line (starts with single '#')
if [ "$KEEP_PASS_ONLY" -eq 1 ]; then
	filter_params=$filter_params" | grep -E 'PASS|^#'"
fi
if [ "$FILTER_SEGDUP" -eq 1 ]; then
	filter_params=$filter_params" | grep -Ev 'segdup'"
fi
if [ "$FILTER_LCR" -eq 1 ]; then
	filter_params=$filter_params" | grep -Ev 'lcr'"
fi
echo "Filters specified: "$filter_params



# ==================== Main script ====================
echo "chr${chr} (1/5) Retrieving POS, REF, ALT, QUAL ..."
eval "zcat < $vcf_dir/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz | grep -v '##' | $filter_params | cut -f2,4,5,6 > $tmp_dir/chr${chr}_gnomad_subtable.${population}.txt"


echo "chr${chr} (2/5) Retrieving AC, AF, AN (Allele Count, Frequency, Number) and DP (Depth of informative coverage for each sample) ..."
eval "zcat < $vcf_dir/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz | grep -v '##' | $filter_params | cut -f8 | cut -d';' -f1 | sed 's/INFO/AC/' > $tmp_dir/chr${chr}_gnomad_ac.${population}.txt"
eval "zcat < $vcf_dir/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz | grep -v '##' | $filter_params | cut -f8 | cut -d';' -f2 | sed 's/INFO/AN/' > $tmp_dir/chr${chr}_gnomad_an.${population}.txt"
eval "zcat < $vcf_dir/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz | grep -v '##' | $filter_params | cut -f8 | cut -d';' -f3 | sed 's/INFO/AF/' > $tmp_dir/chr${chr}_gnomad_af.${population}.txt"
eval "zcat < $vcf_dir/gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz | grep -v '##' | $filter_params | cut -f8 | sed 's/.*DP=/DP=/' | sed 's/;.*//' | sed 's/INFO/DP/' > $tmp_dir/chr${chr}_gnomad_dp.${population}.txt"


echo "chr${chr} (3/5) Combinining all tables ..."
# Output order of INFO fields: AC, AN, AF, DP
paste $tmp_dir/chr${chr}_gnomad_subtable.${population}.txt $tmp_dir/chr${chr}_gnomad_ac.${population}.txt $tmp_dir/chr${chr}_gnomad_af.${population}.txt $tmp_dir/chr${chr}_gnomad_an.${population}.txt $tmp_dir/chr${chr}_gnomad_dp.${population}.txt > $tmp_dir/chr${chr}_gnomad_table.${population}.txt


echo "chr${chr} (4/5) Expanding multiple variants at same loci ..."
echo "./expand_variant_entries.pl $tmp_dir/chr${chr}_gnomad_table.${population}.txt $out_dir" # output:  $out_dir/chr${chr}_gnomad_table.txt.filtered
./expand_variant_entries.pl $tmp_dir/chr${chr}_gnomad_table.${population}.txt $out_dir # output:  $out_dir/chr${chr}_gnomad_table.txt.filtered


echo "chr${chr} (5/5) Cleaning up temporary files ..."
rm -f $tmp_dir/chr${chr}_*
