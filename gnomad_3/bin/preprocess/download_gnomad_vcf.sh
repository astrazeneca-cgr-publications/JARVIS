#!/bin/bash
#SBATCH -J downl_gnomad_vcf
#SBATCH -o out.downl_gnomad_vcf.txt
#SBATCH -e err.downl_gnomad_vcf.txt
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00 

out_dir="../../vcf" 
mkdir -p $out_dir; cd $out_dir

version="3.0"

for i in `seq 1 22` X;
do
        echo Downloading chr: $i
        wget --no-check-certificate "https://storage.googleapis.com/gnomad-public/release/${version}/vcf/genomes/gnomad.genomes.r${version}.sites.chr${i}.vcf.bgz" &

done

wait
