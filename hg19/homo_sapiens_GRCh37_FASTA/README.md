```
# download fasta files for each chromosome from Ensembl (human build GRCh37 - hg19)
./download.sh

## [IMPORTANT]: replace all header with >chr1, >chr2, etc...
./fix_headers.sh


# concatenate chromosome fasta files into a single one
./concat.sh


# Create index file to extract sequences from 2bit file
./build_index.sh
