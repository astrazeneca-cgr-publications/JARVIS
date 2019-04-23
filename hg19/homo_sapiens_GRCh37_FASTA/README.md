```
# download fasta files for each chromosome from Ensembl (human build GRCh37 - hg19)
./download.sh

## [IMPORTANT]: replace all header with >chr1, >chr2, etc...

# concatenate chromosome fasta files into a single one
./concat.sh

# unzip concatenated fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.all_chromosome.fa.gz

# create sequence dictionary (.dict) file
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh37.75.dna.all_chromosome.fa

# create .fai index file
samtools faidx Homo_sapiens.GRCh37.75.dna.all_chromosome.fa
```

# Retrieve sequences from indexed human genome
``` 
../../bin/util/faToTwoBit Homo_sapiens.GRCh37.75.dna.all_chromosome.fa hsa37.2bit   

../../bin/util/twoBitInfo hsa37.2bit stdout | sort -k2rn > hsa37.chrom.sizes   

## Test: ../../bin/util/twoBitToFa hsa37.2bit:1:1-10 /dev/stdout 
```
