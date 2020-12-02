# Clinvar-denovodb training set coordinates BED file
tail -n+2 ../../../out/topmed-NEW_ClinVar_pathogenic-denovodb_benign-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/clinvar_feature_tables/full_feature_table.clinvar_denovodb.bed | cut -f1-3 > clinvar_denovodb.bed


# Split into multiple files (max. filesize for DeepSEA server: 2MB)
split -l 60000 --additional-suffix=".bed" clinvar_denovodb.bed "clinvar_denovodb."


# Run on DeepSEA web-server; Download results and unzip


# Concatenate DeepSEA scores per dataset into BED files
```
./create_full_deepsea_bed_file.sh
```
