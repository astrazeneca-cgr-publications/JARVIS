import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter
import pandas as pd
import dask.dataframe as dd
import numpy as np
from stats_util import is_outlier
import sys, os
from scipy.stats import mannwhitneyu
from scipy import interp

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)




class SVexplorer:

	def __init__(self, maf=None, af_mode='rare'):

		self.genomic_classes = ['utr', 'lof', 'promoter', 'copy_gain', 'dup_partial', 'inv_span', 'dup_lof', 'intronic', 'intergenic']
		self.maf = maf
		self.af_mode = af_mode


		# ==== Dirs init ====
		self.sv_bed_dir = 'sv-bed'
		

	
	
	def filter_sv_by_maf(self, sv_df):
	   
		if self.maf is not None:
			if self.af_mode == 'rare':
				sv_df = sv_df.loc[ sv_df['af'].astype(float) <= maf]
				
			elif self.af_mode == 'common':
				sv_df = sv_df.loc[ sv_df['af'].astype(float) >= maf]
	
		return sv_df


	
	def process_scores_per_class(self):
			
		total_len_df = pd.DataFrame()
	
		# ======  Read gnomAD SV BED file per Genomic class  ======
		for genomic_class in self.genomic_classes:
		
			print('\n', self.sv_bed_dir + '/SV.' + genomic_class + '.bed')
			sv_df = pd.read_csv(self.sv_bed_dir + '/SV.' + genomic_class + '.bed', sep='\t', header=None, low_memory=False)

			sv_df.columns = ['chr', 'start', 'end', 'af']
			sv_df['chr'] = sv_df['chr'].str.replace('chr', '')
			print('sv_df:', sv_df.shape)
			#print(sv_df.head())
			#print(sv_df.info())

	  
			# Filter Structural Variants by MAF
			sv_df = self.filter_sv_by_maf(sv_df)
			
			print('sv_df - filtered by AF:', sv_df.shape)
			cur_sv_lengths = (sv_df['end'] - sv_df['start']).tolist()
			cur_sv_lengths = np.log10(cur_sv_lengths)
			
			tmp_df = pd.DataFrame(cur_sv_lengths, columns=[genomic_class], index=list(range(len(cur_sv_lengths))))
			
			if total_len_df.shape[0] > 0:
				total_len_df = pd.concat([total_len_df, tmp_df], axis=1)
			else:
				total_len_df = tmp_df
			print('total_len_df:', total_len_df.shape)
			
			
		# sort genomic classes by median
		total_len_df = total_len_df.reindex( total_len_df.median().sort_values().index, axis=1 )

		fig, ax = plt.subplots(figsize=(12, 6))
	
		#total_len_df.boxplot(grid=False)
		
		sns.boxplot(data=total_len_df, linewidth=1, palette="Set2", width=0.5)
		
		plt.title('SV - Length distribution')
		plt.xlabel('Genomic classes')
		plt.ylabel('nt-length distribution\n(log10 scale)')
		plt.show()

		pdf_filename = self.sv_bed_dir + '/SV_len_districution-across_all_classes.pdf'
		fig.savefig(pdf_filename, bbox_inches='tight')



if __name__ == '__main__':
	
	maf = 0.001 # None
	af_mode = 'rare' # common
	
	
	sb = SVexplorer(maf=maf, af_mode=af_mode)
	sb.process_scores_per_class()
   
