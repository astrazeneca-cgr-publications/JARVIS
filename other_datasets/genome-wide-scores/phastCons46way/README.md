# 1. ClinVar
```
tabix All_chromosomes.phastCons46way.bed.bgz -B ../../clinvar/clinvar.pathogenic.bed > phastCons46way.clinvar_pathogenic.bed &
tabix All_chromosomes.phastCons46way.bed.bgz -B ../../clinvar/clinvar.benign.bed > phastCons46way.clinvar_benign.bed &

cat phastCons46way.clinvar_pathogenic.bed phastCons46way.clinvar_benign.bed > phastCons46way.clinvar.bed
````

# 2. HGMD
```
tabix All_chromosomes.phastCons46way.bed.bgz -B ../../hgmd/hgmd.pathogenic.bed > phastCons46way.hgmd_pathogenic.bed &

cat phastCons46way.hgmd_pathogenic.bed phastCons46way.clinvar_benign.bed > phastCons46way.hgmd.bed
```
