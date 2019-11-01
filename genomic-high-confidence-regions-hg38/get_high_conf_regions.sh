#!/bin/bash

gnomad_dir="gnomad_3"
gnomad_high_cov_file="../${gnomad_dir}/coverage-files/high_cov_bed_files-min_depth20/gnomad.whole_genome.high_cov.min_depth20.bed"

# ====================================================================
speedseq_regions="ceph18.b37.exclude.2014-01-15.liftovered_to_hg38.bed"
ucsc_repeat_masker="UCSC_RepeatMasker.hg38.bed.gz"
ucsc_seg_dups="UCSC_SegmentDups.hg38.bed.gz"
# ====================================================================


tmp_dir="tmp"
mkdir -p $tmp_dir

echo "Concatenate Repeat masker and Segment duplication regions"
cat $ucsc_repeat_masker $ucsc_seg_dups > $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.bed.gz  

echo "Replace 'chr' with '' (empty string)"  
zcat < $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.bed.gz | sed 's/chr//g'| gzip > $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.no_chr.bed.gz

echo "Sort and merge the concatenated file" 
zcat < $tmp_dir/cat.UCSC_RepeatMasker.SegmentDups.no_chr.bed.gz | sort -k1,1 -k2,2n | mergeBed -i stdin | gzip > $tmp_dir/merged.cat.UCSC_RepeatMasker.SegmentDups.bed.gz  



echo "Filter out UCSC Repeats and SegDup regions from gnomad cov file with depth >= 20"
subtractBed -a $gnomad_high_cov_file -b $tmp_dir/merged.cat.UCSC_RepeatMasker.SegmentDups.bed.gz | gzip > $tmp_dir/sufficient_cov_minus_UCSC_repeats_and_segdups.bed.gz

echo "Filter out regions with extremely high depth from SpeedSeq" 
subtractBed -a $tmp_dir/sufficient_cov_minus_UCSC_repeats_and_segdups.bed.gz -b $speedseq_regions | gzip > high_conf_genomic_regions.bed.gz


# END: resulting file is the final set of high confidence genomic regions for use with gwRVIS 
