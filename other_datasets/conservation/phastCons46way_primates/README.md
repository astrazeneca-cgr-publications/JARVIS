# Workflow

```
./download_wigFix_files.sh

sbatch ./convert_wig2bed.sh
```

# Intersect phastCons BED files with high_conf_regions
[sbatch] ./intersect_high_conf_regions.sh


# Intersect full feature table with high_conf_regions-phastCons (keeping only phastCons entries defined in the full feature table), and then annotate the expanded full feature table with the most/least conserved regions in each genomic class (might take extra time but should worth for not adding unnecessary complexity in the program)
- Or do not intersect full feature table at this stage and get most/least conserved nt-positions for each genomic clas (not ideal)



# -------------

# Intersect with BED file with all calcualted gwRVIS windows
- Look into: gwrvis_core/annotate_feat_table_w_mut_exl_genomic_class.sh (adding clinvar annotation)

# Intersect with BED file with intergenic, UTRs, lincRNAs, UCNEs and VISTA regions.


sbatch ./sort_bed_files_by_conservation.sh
