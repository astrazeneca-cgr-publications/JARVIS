# Prepare train, validation and test sets for a given config file and top_ratio (of tolerant/intolerant windows)
```
python prepare_data.py ../config.yaml 0.001 0   # last argument is `random_seqs` (0 or 1)
```

# Train model
```
python train_nn_model.py ../config.yaml 0.001 0
```
