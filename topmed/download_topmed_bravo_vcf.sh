#!/bin/bash -l 
#SBATCH -J download_bravo_vcf  # Set name for current job 
#SBATCH -o out.downl_bravo_vcf  # Set job output log 
#SBATCH -e err.downl_bravo_vcf  # Set job error log 
#SBATCH --cpus-per-task=1         # Request 1 CPU (core) on a single node 
#SBATCH --mem=8000          # Request amount of memory 
#SBATCH -t 48:0:0           # Request 48 hours runtime

curl 'https://bravo.sph.umich.edu/freeze5/hg38/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: _sm_au_c=iVVt1pWSkkZtZ5fS0e; remember_token="djifos@gmail.com|ce8e1dd6062a5367cb12b4c6afbd531769b5e69a2e3a5d5b8e3c29f0804d71366b8bb6787c17494f4dadaa0b5c3f64e851dad27429e9b2f341fc10a9f52f308b"' --compressed > vcf/bravo-dbsnp-all.vcf.gz
