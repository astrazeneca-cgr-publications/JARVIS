import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
import pandas as pd
import numpy as np
from stats_util import is_outlier
import sys, os
from scipy.stats import mannwhitneyu
from scipy import interp
from imblearn.under_sampling import RandomUnderSampler


from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_score



class ScoreBenchmarking:

	def __init__(self, base_score, use_gwrvis, pathogenic_class, benign_class='intergenic', maf=None):
	
		self.base_score = base_score
		self.use_gwrvis = use_gwrvis		
		self.pathogenic_class = pathogenic_class
		self.benign_class = benign_class
		self.maf = maf

		if self.use_gwrvis:
			self.score = 'gwrvis'
		else:
			self.score = self.base_score
		print("Benchmarked score:", self.score)
		
		
		self.figs_dir = 'Figs'
		if not os.path.exists(self.figs_dir):
			os.makedirs(self.figs_dir)



	def read_sv_scores(self, sv_class):

		input_dir = self.base_score + '-sv-bed'
		
		print("Reading " + self.score + " for the '" + sv_class + "' class...")
		df = pd.read_csv(input_dir + '/' + self.base_score + '.SV.' + sv_class + '.bed', sep='\t', header=None, low_memory=False)

		if self.base_score == 'jarvis':
			#df.columns = ['svtype', 'gene', 'AF', 'genomic_class', 'gwrvis', 'jarvis']
			df.columns = ['AF', 'gwrvis', 'jarvis']
		else:

			print(self.base_score, type(self.base_score))
			print(df.columns)

			#tmp_cols = []
			#for c in range(df.shape[1]):
			#	tmp_cols.append('col' + str(c))
			
			#print(tmp_cols)
			#df.columns = tmp_cols #self.base_score

			#df.rename({'col'+str(df.shape[1]-1): self.base_score}, axis=1, inplace=True)

			df.columns = ['AF', self.base_score]
				
		print(df.head())
		print(df.shape)

		return df



	def preprocess_scores_per_class(self, lof_df, interg_df):

		print("\nScores pre-processing...")
		
		if self.maf is not None:
		
			lof_df = lof_df.loc[ lof_df.AF <= maf, :]
			interg_df = interg_df.loc[ interg_df.AF <= maf, :]	
		
		self.lof_scores = lof_df[self.score].copy()
		self.interg_scores = interg_df[self.score].copy()

	


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
	
		self.model = LogisticRegression(solver='lbfgs', max_iter=100, C=1.0)
	
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
			rus = RandomUnderSampler(random_state=0, sampling_strategy=class_ratio)
			self.X, self.y = rus.fit_resample(self.X, self.y)
			print('Balanced sets:', sorted(Counter(self.y).items()))
		
		
	
		cv = StratifiedKFold(n_splits=n_splits, shuffle=True)
		
		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)
		metrics_list = []
		
		fig, ax = plt.subplots(figsize=(10, 10))

	
		fold = 0
		for train, test in cv.split(self.X, self.y):
		
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
					 
					 
			fold += 1
			print("Fold " + str(fold) + " - AUC: " + str(roc_auc))
			
						
			
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
		
		
		pdf_filename = self.figs_dir + '/' + self.score + '.' + pathogenic_class + \
						'.AUC_' + str(self.mean_auc) + '.pdf'
						
		
		fig.savefig(pdf_filename, bbox_inches='tight')
		
		
		print('Mean AUC:', self.mean_auc)
		



	def run(self):

		lof_df = self.read_sv_scores(self.pathogenic_class)
		interg_df = self.read_sv_scores(self.benign_class)

		self.preprocess_scores_per_class(lof_df, interg_df)

		#self.run_mann_whitney()
		
		self.train_with_cv()





if __name__ == '__main__':

	# Input args
	base_score = str(sys.argv[1]) #'cadd' #'jarvis'
	#pathogenic_classes = ['lof', 'dup_lof', 'inv_span', 'dup_partial', 'utr', 'copy_gain', 'promoter', 'intronic']
	pathogenic_classes = ['copy_gain', 'promoter', 'intronic']
	
	maf = 0.001 # None
	
	
	use_gwrvis = False
	if base_score == 'gwrvis':
		use_gwrvis = True
		base_score = 'jarvis'


	for pathogenic_class in pathogenic_classes:
		print("\n\n-- Pathogenic class:", pathogenic_class)
		sb = ScoreBenchmarking(base_score, use_gwrvis, pathogenic_class=pathogenic_class, maf=maf)
		sb.run()


	

	# - Check associations with commont variants that are in LD with certain SVs.
	# We identified 15,634 common SVs (AF > 1%) in strong LD (R2 ≥ 0.8) with at least one common SNV or indel

	# TODO
	# - Look into biorxiv preprint for additional dataset(?)
	# "Examples of SVs with strong noncoding effects are scarce in humans and model organisms,39–41"


	# Run Fisher's exact test for enrhichment of JARVIS sub-windows that have a score > 0.8 vs < 0.2
