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
from scipy import interp
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold, train_test_split, cross_val_score
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params



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

		base_benchmark_dir = self.out_dir + '/clinvar_scores_benchmarking'
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
		cmd = 'intersectBed -a ' + self.out_dir + '/BED/full_genome.All_genomic_classes.bed -b ../other_datasets/' + variant_datasets[variant_type] + '/' + variant_datasets[variant_type] + '.' + variant_type + '.bed | ' + """awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$4"\t"$5}' """ + ' > ' + gwrvis_clinvar_subset_file
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
		
			
		clinvar_subset_file = self.benchmark_dir + '/' + self.score + '.' + variant_datasets[variant_type] + '_' + variant_type + '.bed'
		score_dir = '../other_datasets/genome-wide-scores/' + self.primary_score

		# Keeping only the intervals of the score-of-interest that overlap with gwRVIS annotations (no additional regions from gwRVIS are inserted)
		cmd = 'intersectBed -wo -a ' + score_dir + '/' + self.primary_score + '.' + variant_datasets[variant_type] + '_' + variant_type + '.bed -b ' + gwrvis_clinvar_subset_file + ' | cut -f1,2,3,4,9,10 > ' + clinvar_subset_file
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

		# DEBUG: Uncomment when I'll have more data from HGMD	
		#self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size=0.2, random_state=0)
		
		# For Cross-Validation
		self.X_train = X
		self.y_train = y
		self.X_test = X
		self.y_test = y

		# logistic regression
		model = LogisticRegression(C=1e9, solver='lbfgs', max_iter=10000)
		cv_auc_scores = cross_val_score(model, X, y, cv=5)
		print(cv_auc_scores)
		mean_auc = np.mean(cv_auc_scores)
		print('5-fold CV AUC:', mean_auc)

		model = LogisticRegressionCV(cv=5, random_state=42, solver='lbfgs', max_iter=10000)
		model.fit(self.X_train, self.y_train)

		self.pred_proba = model.predict_proba(self.X_test)[:, 1]

		return model 



	def plot_roc_curve(self, model, genomic_class):


		fpr, tpr, thresholds = roc_curve(self.y_test, self.pred_proba)
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


		
	def depr_run_logistic_regression(self, pathogenic_df, benign_df, genomic_class):

		pathogenic_df['pathogenic'] = 1
		benign_df['pathogenic'] = 0

		df = pd.concat([pathogenic_df, benign_df], axis=0)


		# Logistic Regression
		model = self.fit_logit_regression(df)

		roc_fig, roc_auc = self.plot_roc_curve(model, genomic_class)
		print('==============================')


		return roc_fig, roc_auc

	
	
		

	def run_logistic_regression(self, pathogenic_df, benign_df, genomic_class):

		pathogenic_df['pathogenic'] = 1
		benign_df['pathogenic'] = 0

		df = pd.concat([pathogenic_df, benign_df], axis=0)

		
		if self.multiple:
			print('> Fitting multiple logistic regression...')
			X = df[['score', 'gwrvis']].values
		else:
			X = df['score'].values.reshape(-1, 1)
		y = df['pathogenic'].values
		

		
		
		# Run classifier with cross-validation and plot ROC curves
		cv = StratifiedKFold(n_splits=5)
		classifier = LogisticRegression(C=1e9, solver='lbfgs', max_iter=10000)

		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)

		
		roc_fig, ax = plt.subplots(figsize=(10, 10))
	
		i = 0
		for train, test in cv.split(X, y):
			probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
			# Compute ROC curve and area the curve
			fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
			tprs.append(interp(mean_fpr, fpr, tpr))
			tprs[-1][0] = 0.0
			roc_auc = auc(fpr, tpr)
			aucs.append(roc_auc)
			plt.plot(fpr, tpr, lw=1, alpha=0.3,
					 label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

			i += 1
		plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
				 label='Chance', alpha=.8)

		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)
		plt.plot(mean_fpr, mean_tpr, color='b',
				 label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
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
		plt.title('Receiver operating characteristic example')
		plt.legend(loc="lower right")
		plt.show()
		plt.close()
		
		
		print('Mean AUC:', mean_auc)
		self.roc_curve_data_per_class[genomic_class] = [mean_auc, mean_fpr, mean_tpr]
		
		
		
		print('==============================')

		return roc_fig, mean_auc

		
		

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

		print('===================================')
		print('pathogenic:', set(pathogenic_df['genomic_class']))
		print('benign:', set(benign_df['genomic_class']))		
		print('===================================')

		genomic_classes = list(set(pathogenic_df['genomic_class']) & set(benign_df['genomic_class']))
		genomic_classes = [g for g in genomic_classes if g != 'vista']


		for genomic_class in genomic_classes:
			print('\nGenomic class:', genomic_class)
			try:	
				pathogenic = pathogenic_df.loc[ pathogenic_df['genomic_class'] == genomic_class, ['score', 'gwrvis'] ]
				benign = benign_df.loc[ benign_df['genomic_class'] == genomic_class, ['score', 'gwrvis'] ]
				print('Pathogenic:', len(pathogenic))
				print('Benign:', len(benign))
			except Exception as e:
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

	all_bennchmark_dir = out_dir + '/clinvar_scores_benchmarking/all-scores'
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




	# Plot ROC curves and Density plots across scores per genomic class
	for genomic_class in genomic_classes:
		#fig, ax = plt.subplots(figsize=(10, 10))
		fig = plt.figure(figsize=(12, 8))
		ax1 = plt.subplot2grid((3, 5), (0, 0), colspan=3, rowspan=3)

		dens_row_index = 0
		dens_col_index = 3
		for score, auc in ordered_scores_per_genomic_class[genomic_class]:
			try:
				# ROC CURVE
				roc_auc, fpr, tpr = roc_curve_data_per_score[score][genomic_class]
				#plt.plot(fpr, tpr, color=class_colors[score],
				#	 lw=1.5, label=score +' (AUC = %0.3f)' % roc_auc)

				ax1.margins(1)
				ax1.plot(fpr, tpr, color=class_colors[score],
					 lw=1.5, label=score +' (AUC = %0.3f)' % roc_auc)


				# Density plots
				if 'gwrvis+' in score:
					ax2 = plt.subplot2grid((3, 5), (dens_row_index, dens_col_index), colspan=1, rowspan=1)
					ax2.set_title(score, fontsize=10)
					ax2.set_axis_off()
				else:
					pathogenic, benign = dens_plot_data_per_score[score][genomic_class]

					ax2 = plt.subplot2grid((3, 5), (dens_row_index, dens_col_index), colspan=1, rowspan=1)
					#fig, ax = plt.subplots(figsize=(10, 10))
					sns.distplot(pathogenic, hist=False, kde=True, label='pathogenic (' + str(len(pathogenic)) + ')', ax=ax2)
					sns.distplot(benign, hist=False, kde=True, label='benign (' + str(len(benign)) + ')', ax=ax2)
					ax2.set_title(score, fontsize=10)
					ax2.legend(fontsize=7)
					ax2.tick_params(axis='both', which='major', labelsize=6) 
					ax2.tick_params(axis='both', which='minor', labelsize=6)
					#plt.close()

				dens_col_index += 1
				if dens_col_index > 4:
					dens_row_index += 1
					dens_col_index = 3

				if dens_row_index > 2:
					dens_row_index = 0

			except Exception as e:
				print('[Warning] Omitting ROC curve for ' + score + 'in genomic class: ' + genomic_class)
				

		ax1.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
		ax1.set_xlim([0.0, 1.0])
		ax1.set_ylim([0.0, 1.05])
		ax1.set_xlabel('False Positive Rate')
		ax1.set_ylabel('True Positive Rate')
		ax1.set_title('ROC curves - ' + genomic_class)
		ax1.legend(loc="lower right", fontsize=12)
		#plt.show()
		#plt.close()


		fig.savefig(all_bennchmark_dir + '/' + genomic_class + '.all_scores.pdf', bbox_inches='tight')




if __name__ == '__main__':
	
	config_file = sys.argv[1]


	out_dir = create_out_dir(config_file)
	out_dir = out_dir + '/full_genome_out'
	#out_dir = '../' + out_dir + '/full_genome_out'

	config_params = get_config_params(config_file)
	pathogenic_set = config_params['pathogenic_set']
	benign_set = config_params['benign_set']

	variant_datasets = {'pathogenic': pathogenic_set, 'benign': benign_set}


	print('Pathogenic: ' + pathogenic_set)
	print('Benign: ' + benign_set)

	roc_curve_data_per_score = {}	
	dens_plot_data_per_score = {}

	all_scores = ['phyloP46way', 'phastCons46way', 'orion', 'cadd', 'gwrvis+cadd', 'gwrvis']
	#all_scores = ['gwrvis'] # ['cadd']


	for score in all_scores:
		print('\n\n----------------\n' + score + '\n\n')

		score_obj = ScoreBenchmark(score, out_dir)
		score_obj.run()


		print(score_obj.roc_curve_data_per_class.keys())

		roc_curve_data_per_score[score] = score_obj.roc_curve_data_per_class
		dens_plot_data_per_score[score] = score_obj.dens_plot_data_per_class

		
	# Aggregate results from multiple scores
	plot_multiple_roc_curves(roc_curve_data_per_score, out_dir)
