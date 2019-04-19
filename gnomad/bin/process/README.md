# Workflow  
Overview:
1. First step: `process_all_chr.sh`
	- Output: filtered_variant_tables-[POPULATION_str]-[FILTERS_str] e.g. filtered_variant_tables-all-PASS_ONLY (or similar based on population and FILTER parameters)
	- This script calls get_chr_table.sh for each chromosome
	- It then calls expand_variant_entries.pl to split multiallelic sites into separate lines (keeping the fields AC, AF, AN and DP from the INFO column)

2. Then run: `keep_gwrvis_high_confidence_regions.sh`
	- SCP job submission: `sbatch keep_gwrvis_high_confidence_regions.sh gnomad 1`



## 1. Get filtered VCF tables  
```
./process_all_chr.sh [input_dir_with_vcf_files] [KEEP_PASS_ONLY: 0|1] [FILTER_SEGDUP: 0|1] [FILTER_LCR: 0|1] \
		     [population (optional): leave blank to consider all populations (otherwise: e.g. FIN, ASJ)]
```
#### Examples
```
- All populations and no filters
[sbatch] ./process_all_chr.sh ../../vcf 0 0 0

- Finnish population (and no filters)
[sbatch] ./process_all_chr.sh ../../vcf 0 0 0  FIN
- Askenzi population (and no filters)
[sbatch] ./process_all_chr.sh ../../vcf 0 0 0  ASJ
```

### Complementary population
- To get all data _EXCEPT_ those also existent in a specified population
```
[sbatch] ./complement_process_all_chr.sh vcf 1 FIN     # --> Ouput: vcf-NON_FIN

```

## Getting started
- Run `keep_gwrvis_high_confidence_regions.sh`:
```
[sbatch] keep_gwrvis_high_confidence_regions.sh [dataset: gnomad|topmed] [FILTERED_FLAG: 1|0)

# e.g.
sbatch keep_gwrvis_high_confidence_regions.sh gnomad 1
```





## 2. Retain those with coverage over a threshold
```
# Define input dir:
# 	to choose between all collapsed variant tables --> collapsed_variant_tables 
# 	or collapsed variant tables after filtering variants based on PASS flag and segdup/lcr annotations --> collapsed_variant_tables-FILTERED
#	or VCFs from a particular population (e.g. FIN, ASJ) --> collapsed_variant_tables-FIN-FILTERED
# [ segdup: segment duplication, lcr: low complexity region]

# Default filtering script:

./keep_gnomad_sufficient_coverage_regions.sh [min_depth (e.g. 20)] [FILTER: 1|0] [population (e.g. all, FIN, ASJ)]

> Run:
./keep_gnomad_sufficient_coverage_regions.sh 20 1 FIN
./keep_gnomad_sufficient_coverage_regions.sh 20 1 ASJ

```
```
# ** For getting all data _EXCEPT_ those also existent in a specified population **
./complement_keep_gnomad_sufficient_coverage_regions.sh 20 1 FIN
```



### Appendix
- (For further local-processing):
Download generated output dir to local installation for further runs, e.g.
```
scp -r [user]@[path]:[path]/gnomad_sufficient_coverage_regions_variant_tables-min_depth20-FIN-FILTERED .
```

