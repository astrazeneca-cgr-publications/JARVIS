import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
import sys
import os
from scipy.stats import mannwhitneyu
import random

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir



def calc_mann_whitney_u(a, b):
	# Mann-Whitney U test
	mannwu_pval = mannwhitneyu(a, b).pvalue
	mannwu_pval_str = 'Mann-Whitney U test: ' + str(mannwu_pval)
	print(mannwu_pval_str)

	return mannwu_pval_str


def plot_gwrvis_distr(input_file):

	df = pd.read_csv(input_file, sep='\t', header=None)
	df.columns = ['gwrvis', 'class']

	df['class'] = df['class'].apply(lambda x: 'CCDS' if x.startswith('CCDS') else 'Non_CCDS')
	print(df.head())
	
	ccds_scores = df.loc[ df['class'] == 'CCDS', 'gwrvis'].tolist()
	non_ccds_scores = df.loc[ df['class'] == 'Non_CCDS', 'gwrvis'].tolist()
	sample_ccds_scores = random.sample(ccds_scores, len(non_ccds_scores))
	
	mannwu_pval_str = calc_mann_whitney_u(ccds_scores, non_ccds_scores)
	sample_mannwu_pval_str = calc_mann_whitney_u(sample_ccds_scores, non_ccds_scores)

	mean_ccds = str(round(np.mean(ccds_scores), 3))
	mean_non_ccds = str(round(np.mean(non_ccds_scores), 3))
	mean_sample_ccds = str(round(np.mean(sample_ccds_scores), 3))


	print("Median ccds:", np.median(ccds_scores))
	print("Median non-ccds:", np.median(non_ccds_scores))
	print("Median sample-ccds:", np.median(sample_ccds_scores))


	fig, ax = plt.subplots(figsize=(10, 10))

	sns.distplot(ccds_scores, hist=False, kde=True, label='CCDS (' + str(len(ccds_scores)) + ') Mean: ' + mean_ccds)
	sns.distplot(sample_ccds_scores, hist=False, kde=True, label='Sample CCDS (' + str(len(sample_ccds_scores)) + ') Mean: ' + mean_sample_ccds)
	sns.distplot(non_ccds_scores, hist=False, kde=True, label='Non CCDS (' + str(len(non_ccds_scores)) + ') Mean: ' + mean_non_ccds) 
	

	plt.xlabel('gwRVIS')
	plt.ylabel('Density')
	plt.title('gwRVIS benchmarking on original Orion dataset\n' + mannwu_pval_str + '\n(Sample) ' + sample_mannwu_pval_str + '\nOriginal Orion Mann-Whitney U p-val: 0.001')
	plt.close()

	fig.savefig(out_dir + '/gwRVIS_vs_original_Orion_benchmark.pdf', bbox_inches='tight')



if __name__ == '__main__':

	config_file = sys.argv[1]

	out_dir = create_out_dir(config_file)         
	out_dir = out_dir + '/full_genome_out'
	#out_dir = '../' + out_dir + '/full_genome_out'

	print('out_dir:', out_dir)

	gwrvis_bed = out_dir + '/BED/full_genome.All_genomic_classes.bed'
	orion_original_ccds_bed = '../other_datasets/genome-wide-scores/orion/plos_one_scores/S1DataFile.bed'

	
	intersect_out_file = out_dir + '/BED/gwrvis_vs_original_orion_benchmark.tsv'

	cmd = 'intersectBed -wo -a ' + gwrvis_bed + ' -b ' + orion_original_ccds_bed + ' | cut -f4,10 > ' + intersect_out_file
	print(cmd)
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)                 
	stdout, stderr = p.communicate()
	p.kill() 

	plot_gwrvis_distr(intersect_out_file)
