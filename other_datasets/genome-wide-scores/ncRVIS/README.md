# 1. Download genomic boundaries file
https://doi.org/10.1371/journal.pgen.1005492.s015  -->  journal.pgen.1005492.s015

mv journal.pgen.1005492.s015 ncrvis_genomic_boundaries.txt

# 2. Table with ncRVIS and RVIS scores:
ncRVIS_scores_table.csv


# 3. Convert boundaries files to BED file with ncRVIS annotation (merged)
python parse_genomic_boundaries.py	# --> All_chromosomes.ncRVIS.bed


# 4. Intersect with ClinVar Pathogenic / Benign variant BED files
## a. ClinVar
```
intersectBed -a All_chromosomes.ncRVIS.bed -b ../../clinvar/clinvar.pathogenic.bed > ncRVIS.clinvar_pathogenic.bed &
intersectBed -a All_chromosomes.ncRVIS.bed -b ../../clinvar/clinvar.benign.bed > ncRVIS.clinvar_benign.bed &

cat ncRVIS.clinvar_pathogenic.bed ncRVIS.clinvar_benign.bed > ncRVIS.clinvar.bed
```


## b. HGMD
```
intersectBed -a All_chromosomes.ncRVIS.bed -b ../../hgmd/hgmd.pathogenic.bed > ncRVIS.hgmd_pathogenic.bed &

cat ncRVIS.hgmd_pathogenic.bed ncRVIS.clinvar_benign.bed > ncRVIS.hgmd.bed
```

