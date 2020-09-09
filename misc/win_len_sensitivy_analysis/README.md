### Copy summary AUC/class-size files from runs with different window lengths

- gnomad
```
for W in 500 1000 2000 3000 4000 5000 7000 10000; do cp ../../out/gnomad-regression_beta-winlen_${W}.MAF_0.001.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions/full_genome_out/mean_auc_per_intol_class.W${W}.MAF0.001.txt gnomad/; done
```

- topmed
```
for W in 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000; do cp ../../out/topmed-production-winlen_${W}.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/full_genome_out/mean_auc_per_intol_class.W${W}.MAF0.001.txt topmed/; done
```
