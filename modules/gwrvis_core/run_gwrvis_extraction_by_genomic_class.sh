#!/bin/bash

# file with run parameters
config_log=$1
input_classes=$2


# get path from directory of current script 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


cnt=0
for i in `seq 1 22`;
do
	echo chr: $i
	python $DIR/get_gwrvis_distr_by_genomic_class.py $i $config_log $input_classes &

	if [ $cnt = 8 ]; then
		wait
		cnt=0
	fi
	cnt=$((cnt + 1))
done


#i=X
#echo chr: $i
#python get_gwrvis_distr_by_genomic_class.py $i $config_log $input_classes &
wait
