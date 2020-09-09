#!/bin/bash
#SBATCH -t 24:0:0

#zcat < coord_CDTS_percentile_N11257all.txt.gz | tail -n+2 | cut -f1,2,3,4 | bgzip > coord_CDTS_percentile_N11257all.txt.bgz

tabix -p bed coord_CDTS_percentile_N11257all.txt.bgz
