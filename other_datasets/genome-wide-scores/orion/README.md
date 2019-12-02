# Run 
sbatch ./prepare_variant_specific_files.sh


``` (Deprecated analysis)
sbatch get_orion_midpoints_w_tabix.sh
```

### Dependencies:
- https://github.com/samtools/tabix



# Subset Orion based on convservation scores
```
sbatch ./prepare_conservation_specific_files.sh [labelset_size] [discard_zero_values] 
```

e.g.
sbatch ./prepare_conservation_specific_files.sh 10000 0
