#!/bin/bash

function build_index {
	chr=$1

	bed_dir="bed"
	bgzip $bed_dir/chr${chr}.phastCons46way.primates.high_conf_regions.bed
	tabix -p bed $bed_dir/chr${chr}.phastCons46way.primates.high_conf_regions.bed.gz

}

for chr in `seq 1 22`; do

	echo "Chr $chr"
	build_index $chr &

done

wait
