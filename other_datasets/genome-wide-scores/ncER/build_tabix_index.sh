#!/bin/bash
#SBATCH --time=24:0:0

zcat < ncER_10bpBins_allChr_coordSorted.txt.gz | bgzip > ncER_10bpBins_allChr_coordSorted.txt.bgz

tabix -p bed ncER_10bpBins_allChr_coordSorted.txt.bgz 
