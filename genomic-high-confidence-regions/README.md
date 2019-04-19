## Workflow to get high confidence genomic regions

- Prepare BED file with repeats and segment duplications
```
# Concatenate Repeat and Segment duplication regions
cat UCSC_RepeatMasker.bed.gz UCSC_SegmentDups.bed.gz > cat.UCSC_RepeatMasker.SegmentDups.bed.gz

# > sed 'chr' to ''

# Sort and merge the concatenated file
zcat < cat.UCSC_RepeatMasker.SegmentDups.bed.gz | sort -k1,1 -k2,2n | mergeBed -i stdin | gzip > merged.cat.UCSC_RepeatMasker.SegmentDups.bed.gz
```

- Filter out UCSC Repeats and SegDup regions from gnomad cov file with depth >= 20 
```
subtractBed -a gnomad.whole_genome.high_cov.min_depth20.bed.gz -b merged.cat.UCSC_RepeatMasker.SegmentDups.bed.gz | gzip > sufficient_cov20_minus_UCSC_repeats_and_segdups.bed.gz
```

- Filter out regions with extremely high depth from `SpeedSeq`
Resulting file is the **final set of high confidence genomic regions for use with gwRVIS**.
```
subtractBed -a sufficient_cov20_minus_UCSC_repeats_and_segdups.bed.gz -b ceph18.b37.exclude.2014-01-15.bed | gzip > gwrvis_high_conf_regions_final_set.bed.gz
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
