##
tabix orion.1001.masked.new.txt.gz -B ../../clinvar/clinvar.pathogenic.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > orion.clinvar_pathogenic.bed &
tabix orion.1001.masked.new.txt.gz -B ../../clinvar/clinvar.benign.no_chr_prefix.bed | sed 's/^/chr/g' | cut -f1,2,3,4 > orion.clinvar_benign.bed &

```
sbatch get_orion_midpoints_w_tabix.sh
```

### Dependencies:
- https://github.com/samtools/tabix
