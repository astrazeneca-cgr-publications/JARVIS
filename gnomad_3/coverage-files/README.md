## Get genome-wide mean depth values 
```
zcat < gnomad.genomes.r3.0.coverage.summary.tsv.bgz | tail -n+2 | cut -f2 > gnomad.r3.0.coverage_means.txt
```

## Get mean
```
awkm gnomad.r3.0.coverage_means.txt
```

## Get 10,20,30,40-percentiles of mean coverage depths
```
# [Not-necessary] - current implementation crashes
sbatch ./submit_get_depth_percentile.sh  
```


## Split full coverage tsv file into separate files per chromosome
```
sbatch ./split_coverage_files_by_chr.sh
```
