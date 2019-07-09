import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
import numpy as np
import subprocess
import operator
import sys
import os
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mannwhitneyu
import random

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir




class DenovodbBenchmark:
	
	def __init__(self, full_df, phenotype, score, out_dir, genomic_class_type='gwrvis_class', denovo_db_type='ssc', pubmed_id=''):

		self.full_df = full_df
		self.phenotype = phenotype
		self.score = score
		self.out_dir = out_dir 
		self.denovo_db_type = denovo_db_type
		self.pubmed_id = pubmed_id

		self.genomic_class_type = genomic_class_type # 'gwrvis_class' or 'FunctionClass'
		
		self.roc_curve_data_per_class = {}
		self.dens_plot_data_per_class = {}
		
		self.multiple = False
		self.primary_score = self.score
		if 'gwrvis+' in self.score:
			self.multiple = True
			self.primary_score = self.score.replace('gwrvis+', '')
		
		self.benchmark_dir = self.get_benchmark_dir_for_score(self.denovo_db_type, self.primary_score)


	def get_benchmark_dir_for_score(self, denovo_db_type, score):

		base_benchmark_dir = self.out_dir + '/denovodb_benchmarking-' + denovo_db_type
		if not os.path.exists(base_benchmark_dir):
			os.makedirs(base_benchmark_dir)
			
		pheno_benchmark_dir = base_benchmark_dir + '/' + self.phenotype
		if not os.path.exists(pheno_benchmark_dir):
			os.makedirs(pheno_benchmark_dir)

		benchmark_dir = pheno_benchmark_dir + '/' + score
		if not os.path.exists(benchmark_dir):
			os.makedirs(benchmark_dir)

		return benchmark_dir
		
		
	def subset_original_df_by_phenotype_and_pubmedId(self):
	
		selected_columns = [self.score, 'gwrvis_class', 'PrimaryPhenotype', 'FunctionClass']
		self.df = self.full_df[selected_columns].copy()
		
		# subset by phenotype
		self.df = self.df.loc[ self.df['PrimaryPhenotype'].isin(['control', self.phenotype]), :]
		
		self.df['disease_phenotype'] = -1
		self.df.loc[ self.df['PrimaryPhenotype'] == 'control', 'disease_phenotype'] = 0
		self.df.loc[ self.df['PrimaryPhenotype'] != 'control', 'disease_phenotype'] = 1
		
		# drop rows with NAs (-1) for the integrated denovo-db scores
		if self.score != 'gwrvis':
			self.df = self.df.loc[ self.df[self.score] != -1, :]
		
		# TODO: add filter by pubmed_id
		pass
		
		
	def plot_score_densities(self, cases_df, controls_df, genomic_class):
	
		cases = cases_df.loc[:, self.score].tolist()
		controls = controls_df.loc[:, self.score].tolist()
		
		dens_fig, ax = plt.subplots(figsize=(10, 10))

		sns.distplot(cases, hist=False, kde=True, label='cases (' + str(len(cases)) + ')')
		sns.distplot(controls, hist=False, kde=True, label='controls (' + str(len(controls)) + ')')
		plt.title(genomic_class)
		plt.close()

		self.dens_plot_data_per_class[genomic_class] = [cases, controls]
		
		return dens_fig

		
	def fit_logit_regression(self, df):
	
		if self.multiple:
			print('> Fitting multiple logistic regression...')
			X = df[[self.score, 'gwrvis']].values
		else:
			X = df[self.score].values.reshape(-1, 1)
		y = df['disease_phenotype'].values

		# logistic regression
		model = LogisticRegression(C=1e9, solver='lbfgs')
		model.fit(X, y)

		return model, X, y



	def plot_roc_curve(self, model, df, X, genomic_class):

		df['pred_gene_class_prob'] = model.predict_proba(X)[:, 1]

		fpr, tpr, thresholds = roc_curve(df['disease_phenotype'], df['pred_gene_class_prob'])
		roc_auc = round(auc(fpr, tpr), 3)
		print("AUC : %f" % roc_auc)

		# Plot ROC curve
		fig, ax = plt.subplots(figsize=(10, 10))
		plt.plot(fpr, tpr, color='darkorange',
			 lw=2, label='ROC curve (area = %0.3f)' % roc_auc)
		plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.05])
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.title('Receiver operating characteristic - ' + genomic_class)
		plt.legend(loc="lower right", fontsize=14)
		plt.show()
		plt.close()

		self.roc_curve_data_per_class[genomic_class] = [roc_auc, fpr, tpr]

		return fig, roc_auc



	def run_logistic_regression(self, cases_df, controls_df, genomic_class):

		df = pd.concat([cases_df, controls_df], axis=0)
		df = df.sample(frac=1)

		# Logistic Regression
		model, X, y = self.fit_logit_regression(df)

		roc_fig, roc_auc = self.plot_roc_curve(model, df, X, genomic_class)

		return roc_fig, roc_auc

		

	def run(self):

		# subset table by phenotype and pubmed_id (if applicabel)
		self.subset_original_df_by_phenotype_and_pubmedId()
		
		genomic_classes = list(self.df[self.genomic_class_type].unique())
		genomic_classes.append('All_classes')
		
		for genomic_class in genomic_classes:
			print('\n--- Genomic class:', genomic_class)
			
			cases_df = self.df.loc[ self.df['disease_phenotype'] == 1, :]
			controls_df = self.df.loc[ self.df['disease_phenotype'] == 0, :]
			
			if genomic_class != 'All_classes':
				cases_df = cases_df.loc[ self.df[self.genomic_class_type] == genomic_class, :]
				controls_df = controls_df.loc[ self.df[self.genomic_class_type] == genomic_class, :]
			
			
			try:
				dens_fig = self.plot_score_densities(cases_df, controls_df, genomic_class)
				roc_fig, roc_auc = self.run_logistic_regression(cases_df, controls_df, genomic_class)
			except:
				print('Insufficient data points for genomic class:', genomic_class)
				continue
				
			pp = PdfPages(self.benchmark_dir + '/' + self.score + '.Logistic_Regression_ROC.AUC_' + str(roc_auc) + '.' + genomic_class + '.pdf')
			pp.savefig(roc_fig)
			pp.savefig(dens_fig)
			pp.close()
			
			
			
			
			
			
			
		
		
	
def plot_multiple_roc_curves(denovo_db_type, roc_curve_data_per_score, out_dir, phenotype):

	all_bennchmark_dir = out_dir + '/denovodb_benchmarking-' + denovo_db_type + '/' + phenotype + '/all-scores'
	if not os.path.exists(all_bennchmark_dir):
		os.makedirs(all_bennchmark_dir)

	gwrvis_color = '#e31a1c'
	colors = [c for c in sns.color_palette("Paired", 12).as_hex() if c != gwrvis_color]	
	
	# get all genomic classes	
	genomic_classes = []
	for score in roc_curve_data_per_score.keys():
		genomic_classes.extend(list(roc_curve_data_per_score[score].keys()))
	genomic_classes = list(set(genomic_classes))


	# get scores ordered by AUC for each genomic class
	ordered_scores_per_genomic_class = {}
	for genomic_class in genomic_classes:
		tmp_dict = {}
		for score in roc_curve_data_per_score.keys():

			try:
				roc_auc, _, _ = roc_curve_data_per_score[score][genomic_class]
				tmp_dict[score] = roc_auc
			except:
				print('[Warning] No ' + genomic_class + ' elements available for ' + score + ' score - Skipped.')

		sorted_scores = sorted(tmp_dict.items(), key=operator.itemgetter(1), reverse=True)
		ordered_scores_per_genomic_class[genomic_class] = sorted_scores


	all_scores = list(roc_curve_data_per_score.keys())
	class_colors = dict( zip(all_scores, colors[:len(all_scores)]) )
	class_colors['gwrvis'] = gwrvis_color

	# Plot ROC curves across scores per genomic class
	for genomic_class in genomic_classes:
		fig, ax = plt.subplots(figsize=(10, 10))

		for score, auc in ordered_scores_per_genomic_class[genomic_class]:
			try:
				roc_auc, fpr, tpr = roc_curve_data_per_score[score][genomic_class]
				
				plt.plot(fpr, tpr, color=class_colors[score],
					 lw=1.5, label=score +' (AUC = %0.3f)' % roc_auc)
			except Exception as e:
				print('[Warning] Omitting ROC curve for ' + score + 'in genomic class: ' + genomic_class)
				

		plt.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.05])
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.title('Receiver operating characteristic - ' + genomic_class)
		plt.legend(loc="lower right", fontsize=14)
		plt.show()
		plt.close()

		fig.savefig(all_bennchmark_dir + '/' + genomic_class + '.all_scores.pdf', bbox_inches='tight')
	
	
	
		
		

def intersect_gwrvis_with_denovodb(denovo_db_type='ssc'):

	denovo_db_files = {'ssc': '../../other_datasets/denovo-db/bed/denovo-db.ssc-samples.variants.bed',
			   'non-ssc': '../../other_datasets/denovo-db/bed/denovo-db.non-ssc-samples.variants.bed'}

	denovo_db_file = denovo_db_files[denovo_db_type]
	print(denovo_db_file)

	intersect_out_file = out_dir + '/BED/gwrvis_vs_' + denovo_db_type + '.merged.bed'

	cmd = 'intersectBed -wo -a ' + gwrvis_bed + ' -b ' + denovo_db_file + ' > ' + intersect_out_file
	print(cmd)
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)                 
	stdout, stderr = p.communicate()
	p.kill() 

	out_df = pd.read_csv(intersect_out_file, sep='\t', header=None)
	out_df.columns = ['chr', 'gwrvis_start', 'gwrvis_end', 'gwrvis', 'gwrvis_class', 'denovodb_chr', 'denovodb_start', 'denovodb_end', 'PubmedID', 'PrimaryPhenotype', 'Gene', 'Variant', 'FunctionClass', 'NumProbands', 'NumControls', 'SequenceType', 'Transcript', 'codingDnaSize', 'PolyPhen(HDiv)', 'PolyPhen(HVar)', 'SiftScore', 'CaddScore', 'LofScore', 'LrtScore', 'aux']
	
	out_df.drop(['denovodb_chr'], axis=1, inplace=True)
	print(out_df.head())
	print(out_df['PubmedID'].unique())

	return(out_df)
	
	
	
def calc_mann_whitney_u(a, b):
	# Mann-Whitney U test
	mannwu_pval = mannwhitneyu(a, b).pvalue
	mannwu_pval_str = 'Mann-Whitney U test: ' + str(mannwu_pval)
	print(mannwu_pval_str)

	return mannwu_pval_str



def	run_benchmark(denovo_db_type, phenotype, genomic_class_type):
	
	df = intersect_gwrvis_with_denovodb(denovo_db_type=denovo_db_type)

	roc_curve_data_per_score = {}	
	dens_plot_data_per_score = {}

	
	all_scores = ['gwrvis', 'CaddScore', 'PolyPhen(HDiv)', 'PolyPhen(HVar)', 'SiftScore', 'LofScore', 'LrtScore']
	
	phenotype = 'autism'
	for score in all_scores:
		print('\n\n> ', score)
		ssc_obj = DenovodbBenchmark(df, phenotype, score, out_dir, genomic_class_type=genomic_class_type, denovo_db_type=denovo_db_type)
		ssc_obj.run()
		
		roc_curve_data_per_score[score] = ssc_obj.roc_curve_data_per_class
		dens_plot_data_per_score[score] = ssc_obj.dens_plot_data_per_class

	print(roc_curve_data_per_score.keys())
	print(dens_plot_data_per_score.keys())
	
	plot_multiple_roc_curves(denovo_db_type, roc_curve_data_per_score, out_dir, phenotype)

	



if __name__ == '__main__':

	config_file = sys.argv[1]

	base_out_dir = create_out_dir(config_file)         
	out_dir = '../' + base_out_dir + '/full_genome_out'
	gwrvis_bed = out_dir + '/BED/full_genome.All_genomic_classes.bed'
	#gwrvis_bed = base_out_dir + '/gwrvis_scores/full_genome.all_gwrvis.no_win_index.bed'
		
	
	phenotype = 'autism'
	genomic_class_type = 'gwrvis_class' # 'FunctionClass' or 'gwrvis_class'
	
	
	# SSC
	denovo_db_type = 'ssc'
	run_benchmark(denovo_db_type, phenotype, genomic_class_type)
	
	
	# Non-SSC
	#denovo_db_type = 'non-ssc'
	#run_benchmark(denovo_db_type, phenotype, genomic_class_type)

	
