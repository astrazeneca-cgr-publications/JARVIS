#!/bin/sh

way=$1 #46 or 100
phastCons_dir=phastCons${way}way
mkdir -p $phastCons_dir
cd $phastCons_dir


URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons${way}way"
if [ $way == 46 ]; then
	URL="$URL/vertebrate"
elif [ $way == 100 ]; then 
	URL="$URL/hg19.100way.phastCons"
fi


cnt=1
#for chr in `seq 1 22`; do
for chr in 3 11; do

	wget $URL/chr${chr}.phastCons${way}way.wigFix.gz &

	if [ $cnt == 7 ]; then
		wait
		cnt=1
	fi
	cnt=$((cnt+1))
	echo $cnt
done

wait
