#!/bin/bash
#SBATCH -J downl_gnomad_vcf
#SBATCH -o out.downl_gnomad_vcf.txt
#SBATCH -e err.downl_gnomad_vcf.txt
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00 

out_dir="../vcf" 
mkdir -p $out_dir; cd $out_dir

for i in `seq 1 22`;
do
        echo Downloading chr: $i
        wget --no-check-certificate "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.${i}.vcf.bgz" &
done


echo Downloading chr: X
wget --no-check-certificate "https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.chrX.vcf.bgz" &

wait
