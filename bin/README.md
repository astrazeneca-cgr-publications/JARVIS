# Workflow
### Master script:
```
./wgs.sh config.yaml input_classes.txt
```

### More analytically:

- Read input parameters
```
config_log=config.log;
input_classes=input_classes.txt;  
```
<br>

- Record features across fixed and tiled genomic windows (e.g. common/all variants, mut. rate, CpG islands, GC content, etc.) 
```
./parse_all_chromosomes.sh $config_log;
```
<br>


- Perform logistic regression (common ~ all variants) to get gwRVIS scores 
```
python run_full_regression.py $config_log;   
```
<br>


- Convert window indexes (0-based) to real genomic coordinates 
```
python convert_window_indexes_to_genomic_coords.py $config_log;
```
<br>

- Get gwRVIS distribution by genomic class 
```
./run_gwrvis_extraction_by_genomic_class.sh $config_log $input_classes;
```
<br>

- Aggregate gwRVIS scores from all chromosomes
```
python aggregate_gwrvis_scores.py $config_log $input_classes;
```
<br>
