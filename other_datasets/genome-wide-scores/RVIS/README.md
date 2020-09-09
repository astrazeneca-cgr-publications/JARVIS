- Downloaded RVIS_ExACv2: novel unpublished RVIS gene score based on ExAC v2 release 2.0 (accessed: March 15th 2017). As of this release we use CCDS release 20 and Ensembl release 87 annotations ("RVIS_Unpublished_ExACv2_March2017.txt")

- Retain raw RVIS scores per gene (OE-ratio-[ExAC v2])
```
tail -n+2 RVIS_Unpublished_ExACv2_March2017.txt | cut -f1,3 | sed 's/ /\t/g' > RVIS_scores_per_gene.tsv
```

- Download BED file with gene names (from Ensembl Biomart - hg19): "ccds.gene_names-genomic_coords.gz"

- Process it into a BED file:
zcat < ccds.gene_names-genomic_coords.gz | tail -n+2 | awk '{print $4"\t"$2"\t"$3"\t"$1}' > CCDS.bed

- Overlap with the RVIS_scores_per_gene.tsv and create BED file for RVIS scores (output: "All_chromosomes.RVIS.bed"):
```
python create_bed_for_rvis.py
```


