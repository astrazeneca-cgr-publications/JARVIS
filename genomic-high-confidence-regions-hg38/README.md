## Workflow to get high confidence genomic regions

- Download from UCSC BED files with repeat-masker and segment-duplication regions (hg38)

- Run
```
./get_high_conf_regions.sh
```


----
#### Total final coverage
```
zcat < gwrvis_high_conf_regions_final_set.bed.gz | awk '{s+=$3-$2+1} END {print s}' 
# 1328825246

# gwrvis vs Orion
$ subtractBed -a orion_full_genomic_real_estate.bed.gz -b gwrvis_high_conf_regions_final_set.bed.gz | awk '{s+=$3-$2} END {print s}' 
291375090
$ subtractBed -b orion_full_genomic_real_estate.bed.gz -a gwrvis_high_conf_regions_final_set.bed.gz | awk '{s+=$3-$2} END {print s}' 
235157911
```
