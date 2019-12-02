# Workflow

```
./download_wigFix_files.sh

sbatch ./convert_wig2bed.sh
```

# Intersect phastCons BED files with high_conf_regions
[sbatch] ./intersect_high_conf_regions.sh


./subset_by_genomic_class.sh



# Get subset of top N most-conserved and N least-conserved regions
./label_region_by_conservation.sh 10000 0     (1st option: subset size; 2nd option: discard regions with 0 conservation score)
