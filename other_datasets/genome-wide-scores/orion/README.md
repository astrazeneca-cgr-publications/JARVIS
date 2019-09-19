# 1. ClinVar
```
tabix orion.1001.masked.new.txt.gz -B ../../clinvar/clinvar.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > orion.clinvar_pathogenic.bed &
tabix orion.1001.masked.new.txt.gz -B ../../clinvar/clinvar.benign.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > orion.clinvar_benign.bed &

cat orion.clinvar_pathogenic.bed orion.clinvar_benign.bed > orion.clinvar.bed
```


# 2. HGMD
```
tabix orion.1001.masked.new.txt.gz -B ../../hgmd/hgmd.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > orion.hgmd_pathogenic.bed &

cat orion.hgmd_pathogenic.bed orion.clinvar_benign.bed > orion.hgmd.bed
```





``` (Deprecated)
sbatch get_orion_midpoints_w_tabix.sh
```

### Dependencies:
- https://github.com/samtools/tabix
