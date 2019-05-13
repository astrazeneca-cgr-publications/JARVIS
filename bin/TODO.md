# TODO LIST:
> wgs.sh:
	```
	./parse_all_chromosomes.sh $config_log;   
	python run_full_regression.py $config_log;  
	python convert_window_indexes_to_genomic_coords.py $config_log;
	./run_gwrvis_extraction_by_genomic_class.sh $config_log $input_classes;
	python aggregate_gwrvis_scores.py $config_log $input_classes;
	```

	# ======== DONE up-to-here ========== ^^ 
	- Optional (for inter-classes comparison): Integrate the next 2-3 scripts called below.
	- Then write a script that selects the most intolerant/tolerant windows across the entire genome (high-conf-regions with variant data)


> Run logistic regression between intergenic and ucnes (or other classes) to see if gwRVIS can distinguish them efficiently.

> Try just a few more versions of the logistic regression for gwRVIS calculations, using bins and/or other features (how about negative binomial regression?)

> Add annotation for CpG islands using pre-built Ensembl tracks for hg19.

> Add methylation data for each window. Think about other data types too.

> Add CNNs and RNNs with Keras functional API


-----------------
# DONE:
- Make sure mut_rate, cpg etc. are calculated for the very last window too.
