#!/bin/bash

gnomad_high_cov_file="../gnomad/coverage-files/high_cov_bed_files-min_depth20/gnomad.whole_genome.high_cov.min_depth20.bed"

tmp_dir="tmp"
mkdir -p $tmp_dir

# Concatenate Repeat masker and Segment duplication regions
cat UCSC_RepeatMasker.bed.gz UCSC_SegmentDups.bed.gz > $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.bed.gz  

# Replace 'chr' with '' (empty string)  
zcat < $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.bed.gz | sed 's/chr//g'| gzip > $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.no_chr.bed.gz

# Sort and merge the concatenated file 
zcat < $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.no_chr.bed.gz | sort -k1,1 -k2,2n | mergeBed -i stdin | gzip > $tmp_dir/merged.cat.UCSC_RepeatMasker.SegmentDups.bed.gz  

# Filter out UCSC Repeats and SegDup regions from gnomad cov file with depth >= 20
subtractBed -a $gnomad_high_cov_file -b $tmp_dir/merged.cat.UCSC_RepeatMasker.SegmentDups.bed.gz | gzip > $tmp_dir/sufficient_cov_minus_UCSC_repeats_and_segdups.bed.gz

# Filter out regions with extremely high depth from `SpeedSeq` 
subtractBed -a $tmp_dir/sufficient_cov_minus_UCSC_repeats_and_segdups.bed.gz -b ceph18.b37.exclude.2014-01-15.bed | gzip > high_conf_genomic_regions.bed.gz

# END: resulting file is the final set of high confidence genomic regions for use with gwRVIS 
