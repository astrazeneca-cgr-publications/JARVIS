# Workflow

## Calculate single-nt resolution gwRVIS
```
config_file="conf/topmed/config.single_nt_gwrvis.yaml"
./submit_single_nt_gwrvis.sh $config_file

sbatch ./compile_gwrvis_single_nt_scores.sh $config_file

sbatch ./create_tabix_indexes.sh $config_file
```


### Plot distributions of single-nt gwRVIS across CCDS or all genomic classes
```
cd modules/ad-hoc/
```



## Create single-nt feature table (for calculating the genome-wide single-nt jarvis):
```
sbatch ./gwrvis_core/for_prediction_annotate_feat_table_w_mut_exl_genomic_class.sh conf/topmed/config.NEW_genome_wide_scores.yaml static_files/input_classes.txt
```


## Calculate single-nt JARVIS
```
for i in `seq 1 22`; do sbatch -o chr${i}.jarvis ./jarvis/deep_learn_raw_seq/submit_prediction_per_chr.sh conf/topmed/config.NEW_genome_wide_scores.yaml both $i; done
```


----


## JARVIS Training 
```
# example config file: config.NEW_ClinVar_pathogenic.yaml

python jarvis/deep_learn_raw_seq/prepare_data.py conf/topmed/config.NEW_ClinVar_pathogenic.yaml

./submit_all_jarvis_jobs.sh conf/topmed/config.NEW_ClinVar_pathogenic.yaml 1
```



## JARVIS prediction on test set
```
# example config file: config.NEW_ncER-GWAS-testing.yaml
python jarvis/deep_learn_raw_seq/prepare_data.py conf/topmed/config.NEW_ncER-GWAS-testing.yaml

python -u jarvis/deep_learn_raw_seq/test_jarvis_model.py conf/topmed/config.NEW_ncER-GWAS-testing.yaml both intergenic,utr,lincrna,ucne,vista 0 1
```

 





## Feature table for Machine Learning / Deep Learning
This is compiled in `gwrvis_core/full_table_intersect_regulatory_features.sh`

The final feature table (with genomic classes, clinvar annotation and all associated features) is saved into:
**clinvar_feature_tables/full_feature_table.[pathogenic_set]_[benign_set].bed**

- **phastCons_primate** needs to be imputed with the median in:
	- jarvis/variant_classification/run_variant_classification.py
	- jarvis/deep_learn_raw_seq/prepare_data.py



### Master script:
#### Main config file:
```
conf/topmed/config.NEW_ClinVar_pathogenic.yaml
```


```
# all
[sbatch] ./wgs.sh config.yaml input_classes.txt


# coding
[sbatch] ./wgs.sh config.coding.yaml input_classes.coding.txt

```


#### Run `wgs.sh` for multiple values of MAF, win_len and variant type (SNVs, INDELs or both):
```
./submit_all_jobs.sh
```


### More analytically:

- Read input parameters:
```
config_log=config.yaml;
input_classes=input_classes.txt;  
```
<br>

- Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.):
```
./parse_all_chromosomes.sh $config_log;
```
<br>


- Perform logistic regression (common ~ all variants) to get gwRVIS scores: 
```
python run_full_regression.py $config_log;   
```
<br>


- Convert window indexes (0-based) to real genomic coordinates:
```
python convert_window_indexes_to_genomic_coords.py $config_log;
```
<br>

- Get gwRVIS distribution by genomic class:
```
./run_gwrvis_extraction_by_genomic_class.sh $config_log $input_classes;
```
<br>

- Aggregate gwRVIS scores from all chromosomes:
```
python aggregate_gwrvis_scores.py $config_log $input_classes;
```
<br>


- Get gwRVIS distribution per genomic class across the entire genome:
```
python get_whole_genome_rvis_distr.py $config_log $input_classes;
```
<br>

- Make refined density plots and boxplots for gwRVIS distribution per genomic class with `ggrdiges`:
```
python make_ggridges_plots.py -c $config_log;
```
