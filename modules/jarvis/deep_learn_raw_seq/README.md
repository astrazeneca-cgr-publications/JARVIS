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
