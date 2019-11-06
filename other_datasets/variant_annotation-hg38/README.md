-- variant_annotation-original contains all files as used till 22/10/2019
[Deprecated: - variant_annotation now subtracts clinvar benign variants from each set of pathogenic]


# Create entries for pathogenic/benign (if applicable) variant BED files for each resource, 
# with and without the 'chr' prefix, appending a column with 'Pathogenic' or 'Benign' respectively in the 4th column:

```
e.g.

# HGMD has only annotation for pathogenic variants
hgmd.pathogenic.no_chr_prefix.bed
hgmd.pathogenic.bed

# ClinVar has annotation for both pathogenic and benign
clinvar.pathogenic.bed
clinvar.benign.bed
```
