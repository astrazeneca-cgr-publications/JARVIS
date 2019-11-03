import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import interp
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, auc
from collections import Counter
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir, get_config_params



def is_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
                points = points[:,None]
       
    median = np.median(points, axis=0)   
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    
    med_abs_deviation = np.median(diff)
    if len(points) == 1 and med_abs_deviation == 0.0:
        return np.array([False])
            
    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh



def get_gwrvis_for_genomic_class(genomic_class, class_annot):

	tmp_df = full_gwrvis_df.loc[full_gwrvis_df.loc[:, 'genomic_class'] == genomic_class].copy()

	if discard_positive:
		tmp_df = tmp_df.loc[ tmp_df['gwrvis'] <= 0, ]

	tmp_df['gene_annot'] = class_annot
	del tmp_df['genomic_class']

	return tmp_df


def compile_gwrvis_differential_df(intol_class, toler_class): 
	intolerant_df = get_gwrvis_for_genomic_class(intol_class, 'intolerant')
	tolerant_df = get_gwrvis_for_genomic_class(toler_class, 'tolerant') 

	df = pd.concat([tolerant_df, intolerant_df], axis=0)
	df['idx'] = df.index
	df['gene_class'] = df.gene_annot.map({'tolerant':0, 'intolerant':1})

	return df


	
def run_logit_regression_with_cv(df, intol_class):


	feature_cols = ['gwrvis']
	X = df[feature_cols].values
	y = df['gene_class'].values
	
	#print(X)
	#print(y)
	#print(X.shape)
	#print(y.shape)
	print('class elements:', Counter(y))
	intol_class_elems = dict(Counter(y))[1]
	
	# Run classifier with cross-validation and plot ROC curves
	cv = StratifiedKFold(n_splits=10)
	classifier = LogisticRegression(solver='lbfgs')

	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)

	fig, ax = plt.subplots(figsize=(10, 10))

	i = 0
	for train, test in cv.split(X, y):
		print('> Fold:', str(i+1))
	
		probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)
		plt.plot(fpr, tpr, lw=1, alpha=0.3,
				 label='ROC fold %d (AUC = %0.2f)' % (i+1, roc_auc))

		i += 1
		
	plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

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
	plt.title('ROC curves for gwRVIS predictive power (UCNEs vs Intergenic)\nwith stratified 10-fold cross-validation')
	plt.legend(loc="lower right")
	plt.show()
	
	
	discard_positive_str = ''
	if discard_positive:
		discard_positive_str = '.discard_positive'

	cv_mean_auc = round(np.mean(aucs), 3)
		
	pdf_filename = full_genome_out_dir + '/LR_ROC_gwRVIS.intergenic_vs_' + intol_class + '.AUC_' + str(cv_mean_auc) + discard_positive_str + '.pdf'
	fig.savefig(pdf_filename, bbox_inches='tight')
	
	return cv_mean_auc, intol_class_elems
	
	
	
def run_logit_regression(df):
	"""
		Deprecated
	"""
	print(df.head())
	feature_cols = ['gwrvis']
	X = df[feature_cols]
	y = df['gene_class']
	
	#print(X)
	#print(y)

	# logistic regression
	model = LogisticRegression(solver='lbfgs')
	
	#mean_roc_auc = cross_val_score(model, X, y, cv=10, scoring='roc_auc').mean()
	#print('10-fold CV mean ROC_AUC:', mean_roc_auc)
		
	model.fit(X, y)
		
	return model, X, y



def plot_roc_curve(model, df, X):
	"""
		Deprecated
	"""
	df['pred_gene_class_prob'] = model.predict_proba(X)[:, 1]

	fpr, tpr, thresholds = roc_curve(df['gene_class'], df['pred_gene_class_prob'])
	roc_auc = round(auc(fpr, tpr), 6)
	print("Area under the ROC curve: %f" % roc_auc)

	# Plot ROC curve
	fig, ax = plt.subplots(figsize=(10, 10))
	plt.plot(fpr, tpr, color='darkorange',
		 lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right", fontsize=14)
	plt.show()

	discard_positive_str = ''
	if discard_positive:
		discard_positive_str = '.discard_positive'

	pdf_filename = full_genome_out_dir + '/Logistic_Regression_ROC_gwRVIS.AUC_' + str(roc_auc) + discard_positive_str + '.pdf'

	fig.savefig(pdf_filename, bbox_inches='tight')


	

if __name__ == '__main__':

	config_file = sys.argv[1]	#"../config.yaml"
	discard_positive = int(sys.argv[2])	# 1 or 0: discard or not windows with positive selection (gwRVIS > 0)
	
	out_dir = create_out_dir(config_file)    
	full_genome_out_dir = out_dir + '/full_genome_out'
	#full_genome_out_dir = '../' + out_dir + '/full_genome_out'
	
	run_params = get_config_params(config_file)
	win_len = run_params['win_len']

	full_gwrvis_df = pd.read_csv(full_genome_out_dir + '/Whole_genome_gwRVIS_per_class.csv')
	print(full_gwrvis_df.head())
	print(full_gwrvis_df['genomic_class'].unique())


	#print('Full Df size:', df.shape)
	#df = df.loc[ ~is_outlier(df.loc[:,'gwrvis']), :]
	#print('Df size (no outlier):', df.shape)


	intolerant_classes = [x for x in full_gwrvis_df['genomic_class'].unique() if x != 'intergenic']
	print(intolerant_classes)


	mean_auc_per_intol_class = {}
	
	for intol_class in intolerant_classes:
		print('\nIntolerant class:', intol_class)
		toler_class = 'intergenic'

		if 'config.coding.yaml' in config_file:
			intol_class = 'intolerant'
			toler_class = 'tolerant'

		df = compile_gwrvis_differential_df(intol_class, toler_class)


		# ### Logistic Regression
		#df = df[ ~np.isnan(df.gwrvis) ]
		
		
		cv_mean_auc, intol_class_elems = run_logit_regression_with_cv(df, intol_class)
		
		print("Intergenic vs " + intol_class + ": " + str(cv_mean_auc))
		mean_auc_per_intol_class[intol_class] = [cv_mean_auc, intol_class_elems]

		# [Deprecated]
		#model, X, y = run_logit_regression(df)
		#plot_roc_curve(model, df, X)


	with open(full_genome_out_dir + '/mean_auc_per_intol_class.W' + str(win_len) + '.txt', 'w') as out_fh:
		out_fh.write('genomic_class\tauc\tclass_size\n')
		for k in sorted(mean_auc_per_intol_class.keys()):
			out_fh.write(k + '\t' + str(mean_auc_per_intol_class[k][0]) + '\t' + str(mean_auc_per_intol_class[k][1]) + '\n')
