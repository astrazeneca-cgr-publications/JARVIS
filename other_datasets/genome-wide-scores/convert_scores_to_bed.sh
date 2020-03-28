#!/bin/sh
#SBATCH -t 24:0:0
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4

score=$1	# ${score} or phastCons
way=$2	#46 or 100
input_dir=${score}${way}way

echo "$score ${way}-way"

cnt=1
for chr in `seq 1 22`; do
#for chr in 21; do
	echo "Converting chr $chr ..."
	gunzip -c ${score}${way}way/chr${chr}.${score}${way}way.wigFix.gz | convert2bed -i wig -o bed | cut -f1,2,3,5 > ${score}${way}way/chr${chr}.${score}${way}way.bed &

	if [ $cnt == 4 ]; then
		cnt=1
		wait
	fi
	cnt=$((cnt+1))
done

wait
echo "wigFix to bed conversions complete."


# sCconcatenate all bed files into a single one
cat ${score}${way}way/chr*.${score}${way}way.bed > ${score}${way}way/All_chromosomes.${score}${way}way.bed
