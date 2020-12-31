### Training JARVIS for $genomic_classes using 'structured' features..."
```
use_fixed_cv_batches=0
python -u jarvis/deep_learn_raw_seq/train_nn_model.py $config_file structured $genomic_classes $use_fixed_cv_batches $cv_repeats
```


(Sequences & both features)
--------------------------
### Training JARVIS for $genomic_classes using 'sequences' features..."
```
use_fixed_cv_batches=1
python -u jarvis/deep_learn_raw_seq/train_nn_model.py $config_file sequences $genomic_classes $use_fixed_cv_batches $cv_repeats
```


### Training JARVIS for $genomic_classes using 'both' features..."
```
use_fixed_cv_batches=1
python -u jarvis/deep_learn_raw_seq/train_nn_model.py $config_file both $genomic_classes $use_fixed_cv_batches $cv_repeats
```




## Predict genome wide JARVIS

- Run for each chromosome:
```
e.g. chr 21

sbatch ./jarvis/deep_learn_raw_seq/submit_prediction_per_chr.sh conf/topmed/config.genome_wide_scores_v2.yaml both 21
```

Results are stored in out/[out-folder]/ml_data/jarvis_predictions per chromosome


### Compile all results per chromsome:
```
single-nt-JARVIS/compile_single_nt_jarvis_per_chr.sh
```
