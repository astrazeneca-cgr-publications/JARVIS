#!/bin/bash
#SBATCH -t 24:0:0

mkdir -p eigen-vcf/


echo "Extracting Eigen scores into VCF files from dbNSFP..."
for chr in `seq 1 22`; do

	#tail -n+2 dbNSFPv3.3c/dbNSFP3.3c_variant.chr${chr} | awk -F"\t" '{print $8"\t"$9"\t"$3"\t"$4"\t"$62"\t"$64 }' > eigen-vcf/eigen.chr${chr}.vcf &

	echo $chr
done
wait


echo "Concatenating VCFs from all chromosomes..."
cat eigen-vcf/*.vcf | sort -k1,1V -k2,2n | awk '{if($1 != ".") print $0}' | bgzip > eigen-vcf/eigen.All_chroms.vcf.bgz


echo "Building tabix index for eigen-vcf/eigen.All_chroms.vcf.bgz..."
tabix -p  vcf eigen-vcf/eigen.All_chroms.vcf.bgz
