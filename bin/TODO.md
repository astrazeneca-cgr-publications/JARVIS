# TODO LIST:
> wgs.sh:
	```
	./parse_all_chromosomes.sh $config_log;   
	python run_full_regression.py $config_log;  
	python convert_window_indexes_to_genomic_coords.py $config_log;
	
	```
	# ======== DONE up-to-here ========== ^^ 
	- Integrate the next 3-4 scripts called below.
	- Then write a script that selects the most intolerant/tolerant windows across the entire genome (high-conf-regions with variant data)


> Add annotation for CpG islands using pre-built Ensembl tracks for hg19.

> Add methylation data for each window. Think about other data types too.


-----------------
# DONE:
- Make sure mut_rate, cpg etc. are calculated for the very last window too.
