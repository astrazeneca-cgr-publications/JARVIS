#!/bin/sh
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=22
#SBATCH --time=48:0:0

# -- Keep only coding or non_coding variants from gnomAD into a separate vcf folder

annot=$1 # coding or non_coding

vcf_dir="../../vcf"
ccds_bed="../../../hg19/bed/ccds_from_slave.merged.no_chr_prefix.bed"

vcf_out_dir="../../vcf-${annot}"
mkdir -p $vcf_out_dir
echo "vcf_out_dir: $vcf_out_dir"


function intersect_bed {

	gnomad_chr_file=$1	
	
	header_file=$vcf_out_dir/${gnomad_chr_file}.header
	zcat < $vcf_dir/$gnomad_chr_file | head -2000 | grep '#' > $header_file


	if [ $annot == 'coding' ]; then 
		bedtools intersect -a $vcf_dir/$gnomad_chr_file -b $ccds_bed -wa  > $vcf_out_dir/${gnomad_chr_file}.tmp
	elif [ $annot == 'non_coding' ]; then
		# use -v option with bedtools
		bedtools intersect -v -a $vcf_dir/$gnomad_chr_file -b $ccds_bed -wa  > $vcf_out_dir/${gnomad_chr_file}.tmp
	fi


	cat $header_file $vcf_out_dir/${gnomad_chr_file}.tmp | bgzip -c > $vcf_out_dir/${gnomad_chr_file}
	rm $header_file $vcf_out_dir/${gnomad_chr_file}.tmp 
}


for chr in `seq 1 22`; do
	
	echo $chr
	gnomad_chr_file="gnomad.genomes.r2.1.1.sites.${chr}.vcf.bgz"

	intersect_bed $gnomad_chr_file &

done
wait
