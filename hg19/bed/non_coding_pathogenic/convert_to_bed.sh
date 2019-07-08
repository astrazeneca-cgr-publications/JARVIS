#!/bin/bash

# Convert original Excel entries to BED files
files=`ls *.txt`;
echo $files;

for file in $files;
do	
	echo $file;
	new_name=`echo $file | sed s/\.txt//`;	
	# [TODO: fix deletions] cat $file | awk 'BEGIN{OFS="\t";} NR>1 vv=length($4)-length($3) >= 0? $2+length($4)-length($3) : $2; {print $1,$2,$vv}' > ${new_name}.bed
	cat $file | awk 'BEGIN{OFS="\t";} NR>1 {print $1,$2,$2}' > ${new_name}.bed
done
