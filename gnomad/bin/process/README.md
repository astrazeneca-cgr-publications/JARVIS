# Workflow  
Overview:
1. First step: `process_all_chr.sh`
	- Output: filtered_variant_tables-[POPULATION_str]-[FILTERS_str] e.g. filtered_variant_tables-all-PASS_ONLY (or similar based on population and FILTER parameters)
	- This script calls get_chr_table.sh for each chromosome
	- It then calls expand_variant_entries.pl to split multiallelic sites into separate lines (keeping the fields AC, AF, AN and DP from the INFO column)

2. Then run: `keep_gwrvis_high_conf_regions.sh`
	- SCP job submission: `sbatch keep_gwrvis_high_conf_regions.sh [gnomad|topmed] [input_filtered_dir]`



## 1. Get filtered VCF tables  
```
./process_all_chr.sh [dataset: gnomad|topmed] [vcf_dir: input dir with VCF files] [KEEP_PASS_ONLY: 0|1] [FILTER_SEGDUP: 0|1] [FILTER_LCR: 0|1] \
		     [population (optional): leave blank to consider all populations (otherwise: e.g. FIN, ASJ)]
```
#### Examples
```
- All populations and no filters
[sbatch] ./process_all_chr.sh gnomad ../../vcf 0 0 0

- Finnish population (and no filters)
[sbatch] ./process_all_chr.sh gnomad ../../vcf 0 0 0  FIN
- Askenzi population (and no filters)
[sbatch] ./process_all_chr.sh gnomad ../../vcf 0 0 0  ASJ
```


## 2. Keep filtered variant tables with variants within high confidence genomic regions
- Run `keep_gwrvis_high_conf_regions.sh`:
```
[sbatch] keep_gwrvis_high_conf_regions.sh [dataset: gnomad|topmed] [input_filtered_dir] [population]

e.g.
sbatch keep_gwrvis_high_conf_regions.sh gnomad ../../out/gnomad-filtered_variant_tables-all all
```
