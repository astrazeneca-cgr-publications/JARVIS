#!/bin/bash
#SBATCH -J downl_gnomad_coverage
#SBATCH -o out.downl_gnomad_coverage.txt
#SBATCH -e err.downl_gnomad_coverage.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00

out_dir="../../coverage-files"
mkdir -p $out_dir; cd $out_dir

echo "Downloading gnomAD v3.0 coverage summary TSV ..."
wget --no-check-certificate https://storage.googleapis.com/gnomad-public/release/3.0/coverage/genomes/gnomad.genomes.r3.0.coverage.summary.tsv.bgz
