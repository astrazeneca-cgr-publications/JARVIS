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
tmp_dir = out_dir + '/rvis_distribution/tmp'

tests_dir = out_dir + '/tests_out'
if not os.path.exists(tests_dir):
	os.makedirs(tests_dir)


chroms = list(range(1,23))
print(chroms)


full_rvis_scores = {}

cl = 'phantom-enh'
# depr-temp: cl = 'vista-phantom'
print('>', cl)
full_df = pd.DataFrame()


tmp_out_file = tests_dir + '/full_genome.' + cl + '.csv' 
for chr in chroms:
	#print(chr)	
	tmp_file = tmp_dir + '/rvis_scores_chr' + str(chr) + '.genomic_coords.' + cl + '.bed'
	if not os.path.exists(tmp_file):
		print('[!] No data for this class at chr:', chr)
		continue
	tmp_df = pd.read_table(tmp_file, header=None, sep='\t')
	full_df = pd.concat([full_df, tmp_df], axis=0)

	#print(tmp_df.head())	
	print(tmp_df.shape)	
	#print(full_df.head())
	print(full_df.shape)

full_df.columns = ['chr', 'start', 'end', 'rvis']
full_df.sort_values(by=['rvis'], inplace=True)
full_df.index = range(0, len(full_df))
print(full_df.head())
print(full_df.tail())

full_df.to_csv(tests_dir + '/full_enhancers.bed', header=None, index=False, sep='\t')


## ********** Retrieve top/bottom 20-percentiles from all Phantom5 enhancers **********
top_percentile = 20
top_percentile_ratio = top_percentile / 100.0
lower_bottom_perc_thres = int(top_percentile_ratio * len(full_df))
upper_top_perc_thres = int( (1 - top_percentile_ratio) * len(full_df))

print(lower_bottom_perc_thres, upper_top_perc_thres)

lower_bottom_df = full_df.iloc[ :lower_bottom_perc_thres]
upper_top_df = full_df.iloc[upper_top_perc_thres: ]

print(lower_bottom_df.shape)
print(lower_bottom_df.head())

print(upper_top_df.shape)
print(upper_top_df.head())

lower_bottom_df.to_csv(tests_dir + '/enhancers_lower_bottom_perc.bed', header=None, index=False, sep='\t')
print(tests_dir + '/enhancers_lower_bottom_perc.bed')

upper_top_df.to_csv(tests_dir + '/enhancers_upper_top_perc.bed', header=None, index=False, sep='\t')



fh = open(tests_dir + '/statistical_tests_output.txt', 'w')
# [deprecated] Mann-Whitney U Test between top / lower percentiles
if 0:
	mann_whitn_u_test_str = ''


	a = ubiq_lower_bottom_df['rvis']
	b = ubiq_upper_top_df['rvis']
	#a = lower_bottom_df['rvis']
	#b = upper_top_df['rvis']

	res = mannwhitneyu(a, b)
	print('\n\nUbiquitous Lower ' + str(top_percentile) + '% vs Higher ' + str(top_percentile) + '% of Enhancer RVIS scores:')
	print('>> [MannwhitneyuResult] statistic: ' + str(res.statistic) + ', P-value: ' + str(res.pvalue) + "\r\n")

	fh.write('\n> Ubiquitous Lower ' + str(top_percentile) + '% vs Higher ' + str(top_percentile) + '% of Enhancer RVIS scores:\n')
	fh.write('[MannwhitneyuResult] statistic: ' + str(res.statistic) + ', P-value: ' + str(res.pvalue) + "\r\n")


# get all ubiquitous Phantom5 enhancers
ubiq_phantom5_file = tests_dir + '/ubiquitous_phantom5_enhancers.bed'
cmd = 'intersectBed -a ' + tests_dir + '/full_enhancers.bed  -b ../' + hg_version + '/bed/vista_enhancers_genes_list.bed > ' + ubiq_phantom5_file
p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p1.wait()

# get all non-ubiquitous Phantom5 enhancers
non_ubiq_phantom5_file = tests_dir + '/non_ubiquitous_phantom5_enhancers.bed'
cmd = 'subtractBed -a ' + tests_dir + '/full_enhancers.bed  -b ../' + hg_version + '/bed/vista_enhancers_genes_list.bed > ' + non_ubiq_phantom5_file
p2 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p2.wait()


# perfomr Mann Whitney U test on Phantom5 ubiquitous (VISTA) enhancers vs non-ubiquitous.
ubiq_phantom_df = pd.read_table(ubiq_phantom5_file, header=None)
ubiq_phantom_df.columns = ['chr', 'start', 'end', 'rvis']
print(ubiq_phantom_df.head())

non_ubiq_phantom_df = pd.read_table(non_ubiq_phantom5_file, header=None)
non_ubiq_phantom_df.columns = ['chr', 'start', 'end', 'rvis']
print(non_ubiq_phantom_df.head())

a = ubiq_phantom_df['rvis']
b = non_ubiq_phantom_df['rvis']

print('median ubiquitous:', np.median(a))
print('median non-ubiquitous:', np.median(b))

res = mannwhitneyu(a, b)
print('\n\n Phantom5 ubiquitous vs non-ubiquitous\n')
print('>> [MannwhitneyuResult] statistic: ' + str(res.statistic) + ', P-value: ' + str(res.pvalue) + "\r\n")

fh.write('\n\n Phantom5 ubiquitous vs non-ubiquitous\n')
fh.write('>> [MannwhitneyuResult] statistic: ' + str(res.statistic) + ', P-value: ' + str(res.pvalue) + "\r\n")





# ********** Get ubiquitous lower/upper-top_percentile% enhancers in tests_dir **********
ubiq_lower_bottom_file = tests_dir + '/ubiquitous_lower_bottom_enhancers.bed'
cmd = 'intersectBed -a ' + tests_dir + '/enhancers_lower_bottom_perc.bed  -b ../' + hg_version + '/bed/vista_enhancers_genes_list.bed > ' + ubiq_lower_bottom_file
p1 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p1.wait()

ubiq_upper_top_file = tests_dir + '/ubiquitous_upper_top_enhancers.bed'
cmd = 'intersectBed -a ' + tests_dir + '/enhancers_upper_top_perc.bed  -b ../' + hg_version + '/bed/vista_enhancers_genes_list.bed > ' + ubiq_upper_top_file
p2 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p2.wait()


# Save output log into file
#fh = open(tests_dir + '/statistical_tests_output.txt', 'a')

ubiq_lower_bottom_df = pd.read_csv(ubiq_lower_bottom_file, header=None, sep='\t')
ubiq_lower_bottom_df.columns = ['chr', 'start', 'end', 'rvis']
#print(ubiq_lower_bottom_df.head())
print('Ubiquitous lower-' + str(top_percentile) + '% enhancers:', len(ubiq_lower_bottom_df))
print('All lower-' + str(top_percentile) + '% enhancers:', len(lower_bottom_df))
fh.write('\nUbiquitous lower-' + str(top_percentile) + '% enhancers: ' + str(len(ubiq_lower_bottom_df)) + '\r\n')

ubiq_upper_top_df = pd.read_table(ubiq_upper_top_file, header=None, sep='\t')
ubiq_upper_top_df.columns = ['chr', 'start', 'end', 'rvis']
#print(ubiq_upper_top_df.head())
print('Ubiquitous upper-' + str(top_percentile) + '% enhancers:', len(ubiq_upper_top_df))
print('All upper-' + str(top_percentile) + '% enhancers:', len(upper_top_df))
fh.write('\nUbiquitous upper-' + str(top_percentile) + '% enhancers: ' + str(len(ubiq_upper_top_df)) + '\r\n')

all_lower = len(lower_bottom_df)
ubiq_lower = len(ubiq_lower_bottom_df)
all_upper = len(upper_top_df)
ubiq_upper = len(ubiq_upper_top_df)

contigency_table = [[all_lower, ubiq_lower], [all_upper, ubiq_upper]]
oddsratio, pvalue = stats.fisher_exact(contigency_table)
print("\n>> [Fisher's exact test] odds-ratio: " + str(oddsratio) + ', P-value: ' + str(pvalue))
fh.write("\n[Fisher's exact test] odds-ratio: " + str(oddsratio) + ', P-value: ' + str(pvalue))



#print('Enrichment of ubiquitous enhancers with low RVIS values:', str(float(100 * len(ubiq_lower_bottom_df) / len(ubiq_upper_top_df))) + '%')
perc_enrichment = float(100 * (len(ubiq_lower_bottom_df) - len(ubiq_upper_top_df)) / len(ubiq_upper_top_df))
fold_change = float( len(ubiq_lower_bottom_df) / len(ubiq_upper_top_df) )
print('Enrichment of ubiquitous enhancers with low RVIS values:', str(round(perc_enrichment, 2)) + '%')
fh.write('\nEnrichment of ubiquitous enhancers with low RVIS values: ' + str(round(perc_enrichment, 2)) + '%')

print('Fold-change: x' + str(round(fold_change, 2))) 
fh.write('\nFold-change: x' + str(fold_change)) 
fh.close()
