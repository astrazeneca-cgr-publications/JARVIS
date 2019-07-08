import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import sys, os
import subprocess
import operator
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir



class ScoreBenchmark:
	
	def __init__(self, score, out_dir):

		self.score = score
		self.out_dir = out_dir 

		self.roc_curve_data_per_class = {}
		self.dens_plot_data_per_class = {}
		
		self.multiple = False
		self.primary_score = self.score
		if 'gwrvis+' in self.score:
			self.multiple = True
			self.primary_score = self.score.replace('gwrvis+', '')
		
		self.benchmark_dir = self.get_benchmark_dir_for_score(self.primary_score)



	def get_benchmark_dir_for_score(self, score):

		base_benchmark_dir = self.out_dir + '/scores_benchmarking'
		if not os.path.exists(base_benchmark_dir):
			os.makedirs(base_benchmark_dir)

		benchmark_dir = base_benchmark_dir + '/' + score
		if not os.path.exists(benchmark_dir):
			os.makedirs(benchmark_dir)

		return benchmark_dir



	def intersect_gwrvis_with_clinvar(self, variant_type):
		"""
		    Get gwrivs pathogenic/benign, stratified by genomic class
		"""
		benchmark_dir = self.get_benchmark_dir_for_score('gwrvis')
		
		gwrvis_clinvar_subset_file = benchmark_dir + '/gwrvis.clinvar_' + variant_type + '.bed'
		cmd = 'intersectBed -a ' + self.out_dir + '/BED/full_genome.All_genomic_classes.bed -b ../../other_datasets/clinvar/clinvar.' + variant_type + '.bed | ' + """awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t"$5}' """ + ' > ' + gwrvis_clinvar_subset_file
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr = p.communicate()
		stdout = str(stdout, "utf-8")
		stderr = str(stderr, "utf-8")
		if stderr != '':
			print('[Error]:', stderr)
		p.kill()

		return gwrvis_clinvar_subset_file



	def intersect_clinvar(self, variant_type):
		"""
		    Get other score's pathogenic/benign, stratified by genomic class (intersection with gwrvis results)
		"""
		
		gwrvis_clinvar_subset_file = self.intersect_gwrvis_with_clinvar(variant_type)
		if self.score == 'gwrvis':
			return gwrvis_clinvar_subset_file
		
			
		clinvar_subset_file = self.benchmark_dir + '/' + self.score + '.clinvar_' + variant_type + '.bed'
		score_dir = '../../other_datasets/genome-wide-scores/' + self.primary_score

		# Keeping only the intervals of the score-of-interest that overlap with gwRVIS annotations (no additional regions from gwRVIS are inserted)
		cmd = 'intersectBed -wo -a ' + score_dir + '/' + self.primary_score + '.clinvar_' + variant_type + '.bed -b ' + gwrvis_clinvar_subset_file + ' | cut -f1,2,3,4,9,10 > ' + clinvar_subset_file
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		stdout, stderr = p.communicate()
		stdout = str(stdout, "utf-8")
		stderr = str(stderr, "utf-8")
		if stderr != '':
			print('[Error]:', stderr)
		p.kill()


		return clinvar_subset_file



	def fit_logit_regression(self, df):

		if self.multiple:
			print('> Fitting multiple logistic regression...')
			X = df[['score', 'gwrvis']].values
		else:
			X = df['score'].values.reshape(-1, 1)
		y = df['pathogenic'].values

		# logistic regression
		model = LogisticRegression(C=1e9, solver='lbfgs')
		model.fit(X, y)

		return model, X, y



	def plot_roc_curve(self, model, df, X, genomic_class):

		df['pred_gene_class_prob'] = model.predict_proba(X)[:, 1]

		fpr, tpr, thresholds = roc_curve(df['pathogenic'], df['pred_gene_class_prob'])
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



	def run_logistic_regression(self, pathogenic_df, benign_df, genomic_class):

		pathogenic_df['pathogenic'] = 1
		benign_df['pathogenic'] = 0

		df = pd.concat([pathogenic_df, benign_df], axis=0)

		# Logistic Regression
		model, X, y = self.fit_logit_regression(df)

		roc_fig, roc_auc = self.plot_roc_curve(model, df, X, genomic_class)

		return roc_fig, roc_auc


	def plot_clinvar_densities(self, pathogenic, benign, genomic_class):
		try:
			pathogenic = pathogenic.loc[:, 'score'].tolist()
			benign = benign.loc[:, 'score'].tolist()
			
			fig, ax = plt.subplots(figsize=(10, 10))

			sns.distplot(pathogenic, hist=False, kde=True, label='pathogenic (' + str(len(pathogenic)) + ')')
			sns.distplot(benign, hist=False, kde=True, label='benign (' + str(len(benign)) + ')')
			plt.title(genomic_class)
			plt.close()

			self.dens_plot_data_per_class[genomic_class] = [pathogenic, benign]
			
			return 0, fig

			#fig.savefig(self.benchmark_dir + '/clinvar.pathogenic_vs_benign.density.' + genomic_class + '.pdf')
		except Exception as e:
			print(e)
			print('[Warning]: Insufficient data in genomic class', genomic_class, ' - Skipped.')
			return -1, None


	def run(self):
		pathogenic_file = self.intersect_clinvar('pathogenic')
		benign_file = self.intersect_clinvar('benign')


		pathogenic_df = pd.read_csv(pathogenic_file, header=None, sep='\t')
		benign_df = pd.read_csv(benign_file, header=None, sep='\t')

		pathogenic_df.columns = ['chr', 'start', 'end', 'score', 'gwrvis', 'genomic_class']
		benign_df.columns = ['chr', 'start', 'end', 'score', 'gwrvis', 'genomic_class']

		pathogenic_df.dropna(inplace=True)
		benign_df.dropna(inplace=True)


		genomic_classes = list(set(pathogenic_df['genomic_class']) & set(benign_df['genomic_class']))



		for genomic_class in genomic_classes:
			print('\nGenomic class:', genomic_class)
			try:	
				pathogenic = pathogenic_df.loc[ pathogenic_df['genomic_class'] == genomic_class, ['score', 'gwrvis'] ]
				benign = benign_df.loc[ benign_df['genomic_class'] == genomic_class, ['score', 'gwrvis'] ]
				print('Pathogenic:', len(pathogenic))
				print('Benign:', len(benign))
			except Exception as e:
				print('[Error]:', e)
				print('Insufficient data points for genomic class:', genomic_class)
				continue

			ret, dens_fig = self.plot_clinvar_densities(pathogenic, benign, genomic_class)
			if ret == -1:
				continue

			roc_fig, roc_auc = self.run_logistic_regression(pathogenic, benign, genomic_class)

			pp = PdfPages(self.benchmark_dir + '/' + self.score + '.Logistic_Regression_ROC.AUC_' + str(roc_auc) + '.' + genomic_class + '.pdf')
			pp.savefig(roc_fig)
			pp.savefig(dens_fig)
			pp.close()


def plot_multiple_roc_curves(roc_curve_data_per_score, out_dir):

	all_bennchmark_dir = out_dir + '/scores_benchmarking/all-scores'
	if not os.path.exists(all_bennchmark_dir):
		os.makedirs(all_bennchmark_dir)

	colors = sns.color_palette("Paired", 12).as_hex()
	
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



if __name__ == '__main__':
	
	config_file = sys.argv[1]


	out_dir = create_out_dir(config_file)
	out_dir = '../' + out_dir + '/full_genome_out'

	roc_curve_data_per_score = {}	
	dens_plot_data_per_score = {}

	all_scores = ['phyloP46way', 'phastCons46way', 'orion', 'cadd', 'gwrvis+cadd', 'gwrvis']
	#all_scores = ['orion']


	for score in all_scores:
		print('\n\n----------------\n' + score + '\n\n')

		score_obj = ScoreBenchmark(score, out_dir)
		score_obj.run()


		print(score_obj.roc_curve_data_per_class.keys())

		roc_curve_data_per_score[score] = score_obj.roc_curve_data_per_class
		dens_plot_data_per_score[score] = score_obj.dens_plot_data_per_class

	plot_multiple_roc_curves(roc_curve_data_per_score, out_dir)
