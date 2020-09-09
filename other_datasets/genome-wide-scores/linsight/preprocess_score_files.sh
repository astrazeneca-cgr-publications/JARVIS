#!/bin/bash
#SBATCH -t 24:0:0

#bigWigToWig LINSIGHT.bw linsight.wig
#echo "Successfully converted bigwig to wig"

wig2bed --zero-indexed < linsight.wig | cut -f1,2,3,5 | bgzip > linsight.bed.bgz
echo "Successfully converted wig to bed (zipped with bgzip)"

tabix -p bed linsight.bed.bgz
echo "Successfully built tabix index for LINSIGHT bed"
