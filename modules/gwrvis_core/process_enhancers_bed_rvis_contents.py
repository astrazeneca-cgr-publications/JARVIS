import matplotlib 
matplotlib.use('agg') 
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import pandas as pd
import numpy as np
import sys
from scipy import stats
import subprocess
import os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params


config_file = sys.argv[1]
input_classes_file = sys.argv[2]

run_params = get_config_params(config_file)
hg_version = run_params['hg_version']
grch = {'hg19': '37', 'hg38': '38'}


# Read input classes to aggregate in current call
input_classes = []
with open(input_classes_file) as fh:
	for l in fh:
		if l.startswith('#'):
			continue
		else:
			input_classes.append( l.rstrip() )
print(input_classes)

out_dir = create_out_dir(config_file)
print(out_dir)
tmp_dir = out_dir + '/gwrvis_distribution/BED'

tests_dir = out_dir + '/PHANTOM-tests_out'
if not os.path.exists(tests_dir):
	os.makedirs(tests_dir)


chroms = list(range(1,23))
print(chroms)


full_rvis_scores = {}

cl = 'phantom-enh'
# depr-temp: cl = 'vista-phantom'

def compile_genomic_class_table(cl):
	print('\n> Genomic class:', cl)
	full_df = pd.DataFrame()

	tmp_out_file = tests_dir + '/full_genome.' + cl + '.csv' 
	for chr in chroms:
		#print(chr)	
		tmp_file = tmp_dir + '/gwrvis_scores_chr' + str(chr) + '.genomic_coords.' + cl + '.bed'
		print(tmp_file)

		if not os.path.exists(tmp_file):
			print('[!] No data for this class at chr:', chr)
			continue
		tmp_df = pd.read_csv(tmp_file, header=None, sep='\t')
		full_df = pd.concat([full_df, tmp_df], axis=0)

		#print(tmp_df.head())	
		#print(tmp_df.shape)	
		#print(full_df.head())
		#print(full_df.shape)

	full_df.columns = ['chr', 'start', 'end', 'gwrvis']
	full_df.sort_values(by=['gwrvis'], inplace=True)
	full_df.index = range(0, len(full_df))
	full_df.dropna(inplace=True)

	print(full_df.head())
	print(full_df.tail())
	print(full_df.shape)

	full_df.to_csv(tests_dir + '/' + cl + '.enhancers.bed', header=None, index=False, sep='\t')

	return full_df

phantom_df = compile_genomic_class_table('phantom-enh')
vista_df = compile_genomic_class_table('vista')

vista_std = np.std(vista_df.gwrvis)
vista_mean = np.mean(vista_df.gwrvis)
print("vista_std:", vista_std)
print("vista_mean:", vista_mean)


## ********** Retrieve top/bottom 20-percentiles from all Phantom5 enhancers **********
top_percentile = 10
top_percentile_ratio = top_percentile / 100.0
lower_bottom_perc_thres = int(top_percentile_ratio * len(phantom_df))
upper_top_perc_thres = int( (1 - top_percentile_ratio) * len(phantom_df))

print(lower_bottom_perc_thres, upper_top_perc_thres)

lower_bottom_df = phantom_df.iloc[ :lower_bottom_perc_thres]
upper_top_df = phantom_df.iloc[upper_top_perc_thres: ]

intol_val_upper_thres = lower_bottom_df.tail(1)['gwrvis'].values[0]

print("\n>> Lower bottom:")
print(lower_bottom_df.shape)
#print(lower_bottom_df.head())
#print(lower_bottom_df.tail())
print("intol_val_upper_thres:", intol_val_upper_thres)


print("\n>> Upper top:")
print(upper_top_df.shape)
#print(upper_top_df.head())
#print(upper_top_df.tail())
tol_val_bottom_thres = upper_top_df.head(1)['gwrvis'].values[0]
print("tol_val_bottom_thres:", tol_val_bottom_thres)


lower_bottom_df.to_csv(tests_dir + '/enhancers_lower_bottom_perc.bed', header=None, index=False, sep='\t')
upper_top_df.to_csv(tests_dir + '/enhancers_upper_top_perc.bed', header=None, index=False, sep='\t')





# **************** Get VISTA windows below/above the intolerance/tolerance thresholds to find enrichment through Fisher's exact test ****************

fh = open(tests_dir + '/statistical_tests_output.txt', 'w')


vista_intolerant_df = vista_df.loc[ vista_df.gwrvis < intol_val_upper_thres, :]
vista_tolerant_df = vista_df.loc[ vista_df.gwrvis > tol_val_bottom_thres, :]

print("vista_intolerant_df:", vista_intolerant_df.shape)
print("vista_tolerant_df:", vista_tolerant_df.shape)



all_lower = len(lower_bottom_df)
vista_intolerant = vista_intolerant_df.shape[0]
all_upper = len(upper_top_df)
vista_tolerant = vista_tolerant_df.shape[0]

contigency_table = [[all_lower, vista_intolerant], [all_upper, vista_tolerant]]
oddsratio, pvalue = stats.fisher_exact(contigency_table)
print("\n>> [Fisher's exact test] odds-ratio: " + str(oddsratio) + ', P-value: ' + str(pvalue))
fh.write("\n[Fisher's exact test] odds-ratio: " + str(oddsratio) + ', P-value: ' + str(pvalue))



#print('Enrichment of ubiquitous enhancers with low RVIS values:', str(float(100 * len(ubiq_lower_bottom_df) / len(ubiq_upper_top_df))) + '%')
perc_enrichment = float(100 * (vista_intolerant - vista_tolerant) / vista_tolerant)
fold_change = float( vista_intolerant / vista_tolerant )
print('Enrichment of ubiquitous enhancers with low RVIS values:', str(round(perc_enrichment, 2)) + '%')
fh.write('\nEnrichment of ubiquitous enhancers with low RVIS values: ' + str(round(perc_enrichment, 2)) + '%')

print('Fold-change: x' + str(round(fold_change, 2))) 
fh.write('\nFold-change: x' + str(fold_change)) 
fh.close()
