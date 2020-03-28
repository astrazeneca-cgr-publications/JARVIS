#!/bin/sh

way=$1 #46 or 100
phylop_dir=phyloP${way}way
mkdir -p $phylop_dir
cd $phylop_dir


URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP${way}way"
if [ $way == 46 ]; then
	URL="$URL/vertebrate"
elif [ $way == 100 ]; then 
	URL="$URL/hg19.100way.phyloP100way"
fi


cnt=1
for chr in `seq 1 22`; do

	wget $URL/chr${chr}.phyloP${way}way.wigFix.gz &

	if [ $cnt == 7 ]; then
		wait
		cnt=1
	fi
	cnt=$((cnt+1))
	echo $cnt
done

wait
