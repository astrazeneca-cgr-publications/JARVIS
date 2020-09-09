### [Run]
#### 1.
```
./process_all.sh 1 0  	# (1st arg: keeping PASS only variants; 2nd arg: keeping only SNVs)
```

- Tmp (intermediate) Output: topmed_table.txt
- 15.7% of variants in the entire table do not have a 'PASS' flag in the FILTER field.

- Multiple alleles at the same loci are split in separate lines in the original BRAVO vcf file.


#### 2. 
``` 
[sbatch] keep_gwrvis_high_confidence_regions.sh [input_dir] [dataset; gnomad|topmed] 
# e.g. 
sbatch keep_gwrvis_high_confidence_regions.sh filtered_variant_tables-SNV_only-FILTERED topmed-hg38
```






## Liftover Appendix (hg38 to hg37/19)
- Run these steps before running `process_all.hg19.sh`

#### 1. Download chain files for liftover from:
[http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver)

Required file: `hg38ToHg19.over.chain`

#### 2. Create sequence dictionry for TARGET genome fasta file:
Target genome file: `../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.fa`

```
gatk CreateSequenceDictionary -R ../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.fa

# Output: 
# ../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.dict
```

#### 3. Create samtools fasta index file:
Target genome file: `../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.fa`

```
samtools faidx ../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.fa

# Output: 
# ../hg19/homo_sapiens_GRCh37_FASTA/Homo_sapiens.GRCh37.75.dna.all_chromosome.fa.fai
```
