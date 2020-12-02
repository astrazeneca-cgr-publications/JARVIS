import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
import pandas as pd
import dask.dataframe as dd
import numpy as np
from stats_util import is_outlier
import sys, os
from scipy.stats import mannwhitneyu
from scipy import interp
from multiprocessing import Process, Pool

#from imblearn.under_sampling import RandomUnderSampler

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_score




class ScoreBenchmarking:

	def __init__(self, base_score, genomic_class, maf=None, af_mode='rare'):

		self.base_score = base_score
		self.genomic_class = genomic_class
		self.maf = maf
		self.af_mode = af_mode


		# ==== Dirs init ====
		self.score_bed_dir = self.base_score + '-sv-bed/' + self.genomic_class + '-1_based'
		self.tmp_dir = self.score_bed_dir + '/tmp'
		if not os.path.exists(self.tmp_dir):
			os.makedirs(self.tmp_dir)

		self.processed_scores_file = self.score_bed_dir + "/" + self.base_score + '.SV.' + self.genomic_class + ".processed_scores"
		

	
	
	def filter_sv_by_maf(self, sv_df):
	   
		if self.maf is not None:
			if self.af_mode == 'rare':
				sv_df = sv_df.loc[ sv_df['af'].astype(float) <= maf]
				
			elif self.af_mode == 'common':
				sv_df = sv_df.loc[ sv_df['af'].astype(float) >= maf]
	
		return sv_df


	
	def process_scores_per_interval(self):
		
		# ======  gnomAD SV BED file  ======
		sv_bed_dir = 'sv-bed'
		print(sv_bed_dir + '/SV.' + self.genomic_class + '.bed')
		sv_df = pd.read_csv(sv_bed_dir + '/SV.' + self.genomic_class + '.bed', sep='\t', header=None, low_memory=False)

		sv_df.columns = ['chr', 'start', 'end', 'af']
		sv_df['chr'] = sv_df['chr'].str.replace('chr', '')
		print('sv_df:', sv_df.shape)
		print(sv_df.head())
		print(sv_df.info())

  
		# Filter Structural Variants by MAF
		sv_df = self.filter_sv_by_maf(sv_df)
								


				
				
		# ======  Score Benchmarking  ======
		def get_scores_per_interval(x):

			start = x['start'] + 1
			end = x['end']
  

			elem_scores_df = self.score_df.loc[ start:end ].copy()
			#print(elem_scores_df)
			#print(elem_scores_df.info())
			#print(elem_scores_df.dropna().info())


			elem_scores_df = elem_scores_df.loc[ elem_scores_df['score'] != '.', :]
			elem_scores_df['score'] = pd.to_numeric(elem_scores_df['score'], downcast="float")
			

			mean_score = elem_scores_df['score'].mean(skipna=True)
		
			first_quartile = elem_scores_df['score'].quantile(q=0.25)
			median = elem_scores_df['score'].quantile(q=0.50)
			third_quartile = elem_scores_df['score'].quantile(q=0.75)
			#print('mean:', mean_score)
			

			return pd.Series([mean_score, first_quartile, median, third_quartile])



		def process_scores_per_chr(chrom):

			chrom = str(chrom)

			try:
				self.score_df = pd.read_csv(self.score_bed_dir + '/' + self.base_score + '.SV.' + self.genomic_class + '.with_coords.chr' + str(chrom) + '.1_based', sep='\t', header=None, low_memory=False, index_col=0)
			except:
				print('[Warning]: No data found for chromosome', chrom)
				return

	

			self.score_df.columns = ['score']


			chr_sv_df = sv_df.loc[ sv_df['chr'] == chrom, :]
			print(chr_sv_df.shape)


			chr_mean_scores_df = chr_sv_df.apply(lambda x: get_scores_per_interval(x), axis=1) 
			chr_mean_scores_df.columns = ['mean_score', 'first_quartile', 'median', 'third_quartile']
			
			#print('chr', chrom, ' - chr_mean_scores df:', chr_mean_scores_df.shape)
			chr_mean_scores_df.dropna(inplace=True)
			#print('chr_mean_scores df (without NAs):', chr_mean_scores_df.shape)
			print(chr_mean_scores_df.head())
			print(chr_mean_scores_df.shape)

			
			chr_mean_scores_df.to_csv(self.tmp_dir + '/' + self.base_score + '.SV.' + self.genomic_class + '.chr' + str(chrom) + '.processed_scores', sep='\t', index=False)
			
			print(self.tmp_dir + '/' + self.base_score + '.SV.' + self.genomic_class + '.chr' + str(chrom) + '.processed_scores')


		
		process_jobs = []
		for chrom in range(1,23):
			print('\n>> Chrom:', chrom)
			
			p = Process(target=process_scores_per_chr, args=[chrom])

			process_jobs.append(p)
			p.start()
							
				
		for p in process_jobs:
			p.join()
		


		# Concatenate processed SV scores from all chromosomes
		total_df = pd.DataFrame()
		for chrom in range(1,23):
			cur_file = self.tmp_dir + '/' + self.base_score + '.SV.' + self.genomic_class + '.chr' + str(chrom) + '.processed_scores'
	
			if not os.path.exists(cur_file):
				continue
				
			print(cur_file)
	
			cur_df = pd.read_csv(cur_file, sep='\t')

			if total_df.shape[0] > 0:
				total_df = pd.concat([total_df, cur_df])
			else:
				total_df = cur_df

		total_df.to_csv(self.processed_scores_file, index=False, header=True, sep='\t')
		print('Saved concatenated results to:', self.processed_scores_file)





if __name__ == '__main__':

	base_score = sys.argv[1] #'jarvis'
	genomic_class = sys.argv[2] #'intergenic'
	
	maf = 0.001 # None
	af_mode = 'rare' # common
	
	
	sb = ScoreBenchmarking(base_score, genomic_class=genomic_class, maf=maf, af_mode=af_mode)
	sb.process_scores_per_interval()
   
