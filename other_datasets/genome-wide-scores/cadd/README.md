```
sbatch ./prepare_cadd_raw_scores.sh
```

# Get CADD scores for the subset of Clinvar variants (pathogenic/benign)
```
tabix cadd.whole_genome.all_raw_scores.vcf.bgz -B ../../clinvar/clinvar.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > cadd.clinvar_pathogenic.bed &
tabix cadd.whole_genome.all_raw_scores.vcf.bgz -B ../../clinvar/clinvar.benign.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > cadd.clinvar_benign.bed
```

```
cat cadd.clinvar_pathogenic.bed cadd.clinvar_benign.bed > cadd.clinvar.bed
```
