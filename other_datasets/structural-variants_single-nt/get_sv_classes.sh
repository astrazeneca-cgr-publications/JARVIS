#!/bin/bash

out_dir='sv-bed'
 
zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($6 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.copy_gain.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($7 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.dup_lof.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($8 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.dup_partial.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($10 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.intronic.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($11 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.inv_span.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($12 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.lof.bed &

## (Deprecated - Same as intergenic class)
##zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($14 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.nearest_tss.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($15 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.promoter.bed &

zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if ($4 != "CN=0" && $5 == "True" || $5 == "BOTHSIDES_SUPPORT") print $0}' | awk '{if ($16 != "NA") print $0}' | cut -f1-3,17 | tail -n+2 | sed 's/^/chr/g' > ${out_dir}/SV.utr.bed &
wait



# Create file copies without chr at the beginning of each line
cd $out_dir
sv_classes=('intergenic' 'lof' 'promoter' 'copy_gain' 'utr' 'dup_partial' 'inv_span' 'dup_lof' 'intronic')

for cl in "${sv_classes[@]}"; do
	cat SV.${cl}.bed | sed 's/^chr//g' > SV.${cl}.no_chr.bed
done
