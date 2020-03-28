#!/bin/bash

# Keeping SNPs only
cat denovodb_ctrl_vars.nonSSC.bed | tail -n+2 | awk '{if( length($4) == 1 && length($5) == 1) print $1"\t"$2"\t"$3"\tBenign"}' > denovodb_nonSSC.benign.no_chr_prefix.bed

cat denovodb_ctrl_vars.nonSSC.bed | tail -n+2 | awk '{if( length($4) == 1 && length($5) == 1) print $1"\t"$2"\t"$3"\tBenign"}' | sed 's/^/chr/g' > denovodb_nonSSC.benign.bed
