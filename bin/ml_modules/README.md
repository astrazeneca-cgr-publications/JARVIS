```
cd  [out_dir]/rvis_distribution/tmp

# intolerant;
cat rvis_scores_chr*intolerant.bed > gnomad_intolerant_genes.bed;
cat gnomad_intolerant_genes.bed | grep -v chrX | cut -f4 > gnomad_wgs_intolerant_genes.csv;

# tolerant;
cat rvis_scores_chr*.tolerant.bed > gnomad_tolerant_genes.bed;
cat gnomad_tolerant_genes.bed | grep -v chrX | cut -f4 > gnomad_wgs_tolerant_genes.csv;

# ccds;
cat rvis_scores_chr*.ccds.bed > gnomad_ccds_genes.bed;
cat gnomad_ccds_genes.bed | grep -v chrX | cut -f4 > gnomad_wgs_ccds_genes.csv;

# intergenic;
cat rvis_scores_chr*.intergenic.bed > gnomad_intergenic_genes.bed;
cat gnomad_intergenic_genes.bed | grep -v chrX | cut -f4 > gnomad_wgs_intergenic_genes.csv;

# omim-HI;
cat rvis_scores_chr*.omim-HI.bed > gnomad_omim-HI_genes.bed;
cat gnomad_omim-HI_genes.bed | grep -v chrX | cut -f4 > gnomad_wgs_omim-HI_genes.csv;


cat gnomad_wgs_intolerant_genes.csv gnomad_wgs_tolerant_genes.csv gnomad_wgs_ccds_genes.csv > gnomad.all_ccds_except_omim.csv
```
