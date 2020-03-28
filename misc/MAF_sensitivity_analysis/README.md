# Copy summary AUC/class-size files from runs with different window lengths

W=3000

- gnomad
```
for MAF in 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1; do cp ../../out/gnomad-regression_beta-winlen_${W}.MAF_${MAF}.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions/full_genome_out/mean_auc_per_intol_class.W${W}.MAF${MAF}.txt gnomad/; done
```

- topmed
```
for MAF in 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1; do cp ../../out/topmed-production-winlen_${W}.MAF_${MAF}.varType_snv.Pop_SNV_only-FILTERED/full_genome_out/mean_auc_per_intol_class.W${W}.MAF${MAF}.txt topmed/; done
```
