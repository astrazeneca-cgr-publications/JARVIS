#!/bin/bash
#SBATCH -J downl_gnomad_coverage
#SBATCH -o out.downl_gnomad_coverage.txt
#SBATCH -e err.downl_gnomad_coverage.txt
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00

out_dir="../coverage-files"
mkdir -p $out_dir; cd $out_dir

for i in `seq 1 22`;
do
        echo "Downloading chr: $i ..."
	wget --no-check-certificate "https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/genomes/gnomad.genomes.r2.0.2.chr${i}.coverage.txt.gz" &
done


echo "Downloading chr: X ..."
wget --no-check-certificate "https://storage.googleapis.com/gnomad-public/release/2.0.2/coverage/genomes/gnomad.genomes.r2.0.2.chrX.coverage.txt.gz" &

wait
