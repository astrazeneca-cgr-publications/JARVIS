#!/bin/bash

mkdir -p wigFix

for chr in `seq 1 22`; do

	 wget -P wigFix http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons46way/primates/chr${chr}.phastCons46way.primates.wigFix.gz &
done

wait
