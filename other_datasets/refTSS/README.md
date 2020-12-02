# Get distance from TSS for all ClinVar non-coding pathogenic
```
closestBed -a full_feature_table.clinvar_pathogenic.non_coding_only.sorted.bed -b tss-bed/refTSS_v3.1_human_coordinate.hg19_liftovered.sorted.bed | awk '{ print $1"\t"$2"\t"$3"\t"$4"\t"$2-$6}' > clinvar_pathogenic.distances_from_TSS.csv
```
