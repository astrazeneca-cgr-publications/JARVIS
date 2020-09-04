#!/bin/bash 


config_file=$1 
single_nt_offset=$2   



# ============   Calculate offset limits    ===============
win_len=`cat $config_file | grep "win_len:" | sed 's/^[ \t]*win_len: //' | sed 's/[ ]*$//'`
echo "win_len: **$win_len**"

# single_nt_offset: from -(win_len/2) to +(win_len/2)-1
left_offset=-$((win_len/2))
right_offset=$((win_len/2)) 
right_offset=$((right_offset-1))

# run for one offset only if single_nt_offset arg is set
if [ -n "$single_nt_offset" ]; then
	left_offset="$single_nt_offset"
	right_offset="$single_nt_offset"
fi

echo "left_offset: $left_offset"
echo "right_offset: $right_offset"
# =========================================================


mkdir -p out_logs
mkdir -p err_logs



for single_nt_offset in $( seq $left_offset $right_offset); do 
#for single_nt_offset in $( seq $left_offset 100); do  #DBG

	echo "single_nt_offset: $single_nt_offset"

	out_file="out_logs/single_ntoffset_${single_nt_offset}.out"
	err_file="err_logs/single_ntoffset_${single_nt_offset}.err"

	sbatch -o $out_file -e $err_file ./single_nt_wgs.sh $config_file $single_nt_offset

done


