#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

#import pybedtools
#from imblearn.under_sampling import RandomUnderSampler

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_score


# In[2]:


class ScoreBenchmarking:

	def __init__(self, base_score, pathogenic_class, use_gwrvis=False, benign_class='intergenic', maf=None, af_mode='rare'):

		self.base_score = base_score
		self.use_gwrvis = use_gwrvis		
		self.pathogenic_class = pathogenic_class
		self.benign_class = benign_class
		self.maf = maf
		self.af_mode = af_mode

		if self.use_gwrvis:
			self.score = 'gwrvis'
		else:
			self.score = self.base_score
		print("Benchmarked score:", self.score)


		self.figs_dir = 'Figs'
		if not os.path.exists(self.figs_dir):
			os.makedirs(self.figs_dir)

		# ==== Dirs init ====
		self.score_bed_dir = self.base_score + '-sv-bed/' + pathogenic_class + '-1_based'
		self.tmp_dir = self.score_bed_dir + '/tmp'
		if not os.path.exists(self.tmp_dir):
			os.makedirs(self.tmp_dir)



	def read_sv_scores(self, sv_class):

		input_dir = self.base_score + '-sv-bed'

		print("Reading " + self.score + " for the '" + sv_class + "' class...")
		#df = pd.read_csv(input_dir + '/' + self.base_score + '.SV.' + sv_class + '.bed', sep='\t', header=None, low_memory=False)
		df = dd.read_csv(input_dir + '/' + self.base_score + '.SV.' + sv_class + '.bed', sep='\t', header=None, low_memory=False)



		print(self.base_score, type(self.base_score))
		print(df.columns)

		df.columns = ['AF', self.base_score]

		print(df.head())
		print(df.shape)

		
		return df





	def run_mann_whitney(self):

		print("\nRunning Mann-Whitney U test...")
		res = mannwhitneyu(self.lof_scores, self.interg_scores)
		print('>> [MannwhitneyuResult] ' + self.base_score + ': ' + str(res.statistic) + ', P-value: ' + str(res.pvalue) + "\r\n")




	def plot_densities(self):

		print("\nPlotting densities...")
		fig, ax = plt.subplots(figsize=(12, 12))

		lof_scores = self.lof_scores
		interg_scores = self.interg_scores

		# Remove outliers from gwrvis before plotting
		if self.use_gwrvis:
			lof_scores = self.lof_scores[ ~is_outlier(self.lof_scores, 3.5)]
			interg_scores = self.interg_scores[ ~is_outlier(self.interg_scores, 3.5)]


		lof_scores.plot.kde(bw_method=0.3, color='red')
		interg_scores.plot.kde(bw_method=0.3, color='black')

		fig.savefig(self.figs_dir + '/' + self.score + '_distribution.pdf', bbox_inches='tight')




	def get_metrics(self, test_flat, preds_flat, verbose=0):

		accuracy = accuracy_score(test_flat, preds_flat)
		confus_mat =  confusion_matrix(test_flat, preds_flat)


		TN, FP, FN, TP = confus_mat.ravel()	

		try:
			sensitivity = TP / (TP + FN)
		except:
			sensitivity = 'NA'

		try:	
			precision = TP / (TP + FP)
		except:
			precision = 'NA'

		try:
			specificity = TN / (TN + FP)
		except:
			specificity = 'NA'


		if verbose:
			print('> Confusion matrix:\n', confusion_matrix(test_flat, preds_flat))		
			print('TN:', TN, '\nFP:', FP, '\nFN:', FN, '\nTP:', TP)

			print('\n> Accuracy:', accuracy_score(test_flat, preds_flat))
			print('> Sensitivity:', sensitivity)
			print('> Precision:', precision)
			print('> Specificity:', specificity)

		metrics = {'accuracy': accuracy, 'sensitivity': sensitivity, 'precision': precision, 'specificity': specificity}

		return metrics



	def test_and_evaluate_model(self, y_probas, y_test):


		y_pred_flat = np.argmax(y_probas, axis=1)
		y_test_flat = y_test

		#print('y_test:', y_test_flat)
		#print('y_pred:', y_pred_flat)

		metrics = self.get_metrics(y_test_flat, y_pred_flat)

		return metrics	



	def train_with_cv(self, n_splits=5):

		self.model = LogisticRegression(solver='sag', max_iter=100, C=1.0)
		#self.model = RandomForestClassifier()

		# add pathogenic entries
		X_patho = self.lof_scores.values
		y_patho = np.array([1] * X_patho.shape[0])

		# add control entries (intergenic)
		X_control = self.interg_scores.values
		y_control = np.array([0] * X_control.shape[0])

		self.X = np.append(X_patho, X_control).reshape(-1,1)
		self.y = np.append(y_patho, y_control)


		print(self.X[:10])
		print(self.X.shape)
		print(self.y[:10])
		print(self.y[-10:])
		print(self.y.shape)




		# Fix class imbalance
		positive_set_size = (self.y == 1).sum()
		negative_set_size = (self.y == 0).sum()
		class_ratio = 1/1

		if positive_set_size != negative_set_size:

			print('{ imbalanced sets: ', sorted(Counter(self.y).items()), ' }')
			#rus = RandomUnderSampler(random_state=0, sampling_strategy=class_ratio)
			#self.X, self.y = rus.fit_resample(self.X, self.y)
			
			#X_train, X_test, y_train, y_test = train_test_split( X, y, test_size=0.33, random_state=42, stratify=y)

			
			#print('Balanced sets:', sorted(Counter(self.y).items()))



		cv = StratifiedKFold(n_splits=n_splits, shuffle=True)

		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)
		metrics_list = []

		fig, ax = plt.subplots(figsize=(10, 10))


		fold = 1
		for train, test in cv.split(self.X, self.y):

			print("Fold " + str(fold))
			probas_ = self.model.fit(self.X[train], self.y[train]).predict_proba(self.X[test])			

			# Compute ROC curve and area the curve
			fpr, tpr, thresholds = roc_curve(self.y[test], probas_[:, 1])


			tprs.append(interp(mean_fpr, fpr, tpr))
			tprs[-1][0] = 0.0
			roc_auc = round(auc(fpr, tpr), 3)
			aucs.append(roc_auc)
			plt.plot(fpr, tpr, lw=1, alpha=0.3,
					 label='ROC fold %d (AUC = %0.2f)' % (fold, roc_auc))

			# Evaluate predictions on test and get performance metrics
			metrics = self.test_and_evaluate_model(probas_, self.y[test])
			metrics['auc'] = roc_auc
			metrics_list.append(metrics)		 


			print("Fold " + str(fold) + " - AUC: " + str(roc_auc))
			fold += 1



		plt.plot([0, 1], [0, 1], linestyle='--', lw=1, color='r', label='Chance', alpha=.8)

		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		self.mean_auc = round(auc(mean_fpr, mean_tpr), 3)
		std_auc = np.std(aucs)
		plt.plot(mean_fpr, mean_tpr, color='b',
				 label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (self.mean_auc, std_auc),
				 lw=2, alpha=.8)

		std_tpr = np.std(tprs, axis=0)
		tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
						 label=r'$\pm$ 1 std. dev.')

		plt.xlim([-0.05, 1.05])
		plt.ylim([-0.05, 1.05])
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.title(self.score + ': ' + str(n_splits) + '-fold Cross-Validation ROC Curve')
		plt.legend(loc="lower right")
		plt.show()
		plt.close()


		pdf_filename = self.figs_dir + '/' + self.score + '.' + pathogenic_class +						 '.AUC_' + str(self.mean_auc) + '.pdf'


		fig.savefig(pdf_filename, bbox_inches='tight')


		print('Mean AUC:', self.mean_auc)




	def run(self):

		lof_df = self.read_sv_scores(self.pathogenic_class)
		interg_df = self.read_sv_scores(self.benign_class)

		frac=0.1
		lof_df = lof_df.sample(frac=frac)
		interg_df = interg_df.sample(frac=frac)
		
		self.preprocess_scores_per_class(lof_df, interg_df)

		print('\nComputing lof_scores...')
		
		self.lof_scores = self.lof_scores.compute()
		print(len(self.lof_scores))
		
		
		print('\nComputing interg_scores...')
		self.interg_scores = self.interg_scores.compute()
		
		sample_frac = len(self.lof_scores) / float(len(self.interg_scores))
		self.interg_scores = self.interg_scores.sample(frac=sample_frac)
		print(len(self.interg_scores))
		
		#self.run_mann_whitney()

		self.train_with_cv()
		
	
	
	
	def filter_sv_by_maf(self, sv_df):
	   
		if self.maf is not None:
			if self.af_mode == 'rare':
				sv_df = sv_df.loc[ sv_df['af'].astype(float) <= maf]
				
			elif self.af_mode == 'common':
				sv_df = sv_df.loc[ sv_df['af'].astype(float) >= maf]
	
		return sv_df


	
	def process_scores_per_interval(self, pathogenic_class=None):


		process_scores_file = score_bed_dir + "/" + self.base_score + '.SV.' + pathogenic_class + ".processed_scores", index=False, header=True, sep='\t')
		



		self.cnt = 0
		def get_scores_per_interval(x):
			#print(x)

			start = x['start'] + 1
			end = x['end']
  
			elem_scores_df = self.score_df.loc[ start:end ]		

			mean_score = elem_scores_df['score'].mean(skipna=True)
		
			first_quartile = elem_scores_df['score'].quantile(q=0.25)
			median = elem_scores_df['score'].quantile(q=0.50)
			third_quartile = elem_scores_df['score'].quantile(q=0.75)
			#print('mean:', mean_score)

			return pd.Series([mean_score, first_quartile, median, third_quartile])

		
		
		
		print('\n>> Extracting scores per interval for ==', pathogenic_class, '==')
		if pathogenic_class is None:
			pathogenic_class = self.pathogenic_class
		
		
		# ==== gnomAD SV BED file ====
		sv_bed_dir = 'sv-bed'
		print(sv_bed_dir + '/SV.' + pathogenic_class + '.bed')
		sv_df = pd.read_csv(sv_bed_dir + '/SV.' + pathogenic_class + '.bed', sep='\t', header=None, low_memory=False)

		sv_df.columns = ['chr', 'start', 'end', 'af']
		sv_df['chr'] = sv_df['chr'].str.replace('chr', '')
		print('sv_df:', sv_df.shape)
		print(sv_df.head())
		print(sv_df.info())

  
		# Filter Structural Variants by MAF
		sv_df = self.filter_sv_by_maf(sv_df)
								
				
				

		# ==== Benchmarked score ====
		def process_scores_per_chr(chrom):

			print('\n>> Chrom:', chrom)
			sys.stdout.flush()

			chrom = str(chrom)

			self.score_df = pd.read_csv(score_bed_dir + '/' + self.base_score + '.SV.' + pathogenic_class + 
										'.with_coords.chr' + str(chrom) + '.1_based', 
										sep='\t', header=None, low_memory=False, index_col=0)
			self.score_df.columns = ['score']


			chr_sv_df = sv_df.loc[ sv_df['chr'] == chrom, :]
			#print('\n\n', chr_sv_df.head())


			chr_mean_scores_df = chr_sv_df.apply(lambda x: get_scores_per_interval(x), axis=1) 
			chr_mean_scores_df.columns = ['mean_score', 'first_quartile', 'median', 'third_quartile']
			
			#print('chr', chrom, ' - chr_mean_scores df:', chr_mean_scores_df.shape)
			chr_mean_scores_df.dropna(inplace=True)
			#print('chr_mean_scores df (withouth NAs):', chr_mean_scores_df.shape)
			
			chr_mean_scores_df.to_csv(tmp_dir + '/' + self.base_score + '.SV.' + pathogenic_class + 
										'.chr' + str(chrom) + '.processed_scores', sep='\t', index=False)
			
		"""	
		for chrom in range(1,23):
			
			process_jobs = []
			
			p = Process(target=process_scores_per_chr, args=[chrom])

			process_jobs.append(p)
			p.start()
							
				
		for p in process_jobs:
			p.join()
		"""

		total_df = pd.DataFrame()
		for chrom in range(1,23):
			cur_df = pd.read_csv(tmp_dir + '/' + self.base_score + '.SV.' + pathogenic_class + '.chr' + str(chrom) + '.processed_scores', sep='\t')

			if total_df.shape[0] > 0:
				total_df = pd.concat([total_df, cur_df])
			else:
				total_df = cur_df

		total_df.to_csv(score_bed_dir + "/" + self.base_score + '.SV.' + pathogenic_class + ".processed_scores", index=False, header=True, sep='\t')






if __name__ == '__main__':

	base_score = 'jarvis'
	pathogenic_class = 'intergenic'
	
	maf = 0.001 # None
	af_mode = 'rare' # common
	
	
	sb = ScoreBenchmarking(base_score, pathogenic_class=pathogenic_class, maf=maf, af_mode=af_mode)
	sb.process_scores_per_interval()
   


	# - Check associations with commont variants that are in LD with certain SVs.
	# We identified 15,634 common SVs (AF > 1%) in strong LD (R2 ≥ 0.8) with at least one common SNV or indel

	# TODO
	# - Look into biorxiv preprint for additional dataset(?)
	# "Examples of SVs with strong noncoding effects are scarce in humans and model organisms,39–41"


	# Run Fisher's exact test for enrhichment of JARVIS sub-windows that have a score > 0.8 vs < 0.2


