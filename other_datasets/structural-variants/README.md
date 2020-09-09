# Analaysis steps

#- Get a set of all SVs with a non-NA 'PROTEIN_CODING__LOF' entry (column: 24)
#zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if($9 == "False") print $0}' | awk '{if($12 != "NA" && $4 != "CN=0" && $5 == "True") print $0}' | cut -f1-5,14,17 | sed 's/^/chr/g' > SV.coding_lof.bed


- Get a set of all SVs with 'PROTEIN_CODING__INTERGENIC' = True (column: 21)  
(also excluding one line that has non 'NA' value for 'PROTEIN_CODING__LOF')
```
zcat < gnomad_v2.1_sv.sites.bed.gz | cut -f1-3,5,7,18-28,38 | awk '{if($9 == "True") print $0}' | awk '{if($12 == "NA" && $4 != "CN=0" && $5 == "True") print $0}' | cut -f1-5,14,17 | sed 's/^/chr/g' > SV.intergenic.bed
```

# Get BED files per 'pathogenic'-like SV class (i.e. all except for intergenic)
```
./get_sv_classes.sh
```
