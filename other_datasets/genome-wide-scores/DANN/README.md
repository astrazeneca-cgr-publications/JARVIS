```
sbatch ./prepare_dann_raw_scores.sh
```

# Get DANN scores for the subset of Clinvar variants (pathogenic/benign)
```
tabix DANN_whole_genome_SNVs.tsv.bgz -B ../../clinvar/clinvar.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > dann.clinvar_pathogenic.bed &
tabix DANN_whole_genome_SNVs.tsv.bgz -B ../../clinvar/clinvar.benign.no_chr_prefix.bed | sed 's/^/chr/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$5}' > dann.clinvar_benign.bed
```
