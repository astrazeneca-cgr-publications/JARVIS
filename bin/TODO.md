- wgs.sh:
	```
	./parse_all_chromosomes.sh $config_log;   
	python run_full_regression.py $config_log;  
	```
	# ======== DONE up-to-here ========== ^^ 
	- Run for all chromosomes 
	- Integrate the next 3-4 scripts called below.


- Add methylation data for each window. Think about other data types too.

- Make sure mut_rate, cpg etc. are calculated for the very last window too... [DONE]
