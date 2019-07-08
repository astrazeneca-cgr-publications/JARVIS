All scripts are under `bin/`

## Download gnomAD data

- Download gnomAD VCF data (r2.1.1):
```
[sbatch] ./download_gnomad_vcf.sh     #sbatch for submission to a SLURM-based cluster:

```

- Download gnomAD Coverage data (r2.0.2 -- unchanged in gnomAD r2.1.1):
```
[sbatch] ./download_gnomad_coverage.sh
```


## Calculate simple QC metrics from WGS coverage data
```
run_qc_on_wgs_coverage_data.sh
```

Output:
- `All_mean_depths.gnomad.WGS_coverage.txt` in `../coverage-files`
- `WGS_coverage.QC_depths.txt` with min, max, median and average of mean depths per nt


## Get genomic regions with mean read depth over a threshold
```
[sbatch] ./get_high_cov_gnomad_bed_files.sh [min_depth]
```

Output:
- `high_cov_bed_files-min_depth[min_depth]` under `../coverage-files`



# --------------------------------------------------------
# Filter (retain) coding or non-coding only variants from gnomAD
```
sbatch ./filter_coding_gnomad_regions.sh coding
sbatch ./filter_coding_gnomad_regions.sh non_coding
```
