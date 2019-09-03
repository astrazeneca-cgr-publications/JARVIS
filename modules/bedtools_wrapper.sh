#!/bin/bash

mode=$1
A=$2
B=$3
min_overlap_ratio=$4


# intersect
if [ "$mode" == 'intersect' ]; then
	cat $A | sortBed | bedtools intersect -a stdin -b $B -f $min_overlap_ratio | sortBed | mergeBed -c 4 -o distinct 2>/dev/null
	exit
fi

# intersect -wa [Write the original entry in A for each overlap]
if [ "$mode" == 'intersect_wa' ]; then
	cat $A | sortBed | bedtools intersect -wa -a stdin -b $B -f $min_overlap_ratio | mergeBed -c 4 -o distinct 2>/dev/null
	exit
fi

# subtract
if [ "$mode" == 'subtract' ]; then
	cat $B | sortBed | bedtools subtract -a $A -b stdin -f $min_overlap_ratio | sortBed | mergeBed -c 4 -o distinct
	exit
fi

# subtract -A [Remove entire feature (A) if any overlap with B]
if [ "$mode" == 'subtract_A' ]; then
	cat $B | sortBed | bedtools subtract -A -a $A -b stdin -f $min_overlap_ratio
	exit
fi
