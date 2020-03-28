# Prepare train, validation and test sets for a given config file and top_ratio (of tolerant/intolerant windows)
```
python prepare_data.py ../conf/config.yaml
```

# Train model
```
python train_nn_model.py ../config.yaml 0.001 0
```
