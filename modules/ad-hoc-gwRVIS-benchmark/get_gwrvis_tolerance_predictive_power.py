import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import interp
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, auc
from collections import Counter
import sys
import os
from plot_gwrvis_distr_by_genomic_class import GwrvisDistributionPerClass

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




	
def run_logit_regression_with_cv(df, intol_class):

	df['genomic_class'] = df['genomic_class'].replace({'intergenic': 0, intol_class: 1})
	
	
	feature_cols = ['gwrvis']
	X = df[feature_cols].values
	scaler = StandardScaler()
	#X_std = scaler.fit_transform(X)
	y = df['genomic_class'].values
	
	#print(X)
	#print(y)
	#print(X.shape)
	#print(y.shape)
	print('class elements:', Counter(y))
	intol_class_elems = len(df.loc[ df.genomic_class == 1, :])
	

	# Run classifier with cross-validation and plot ROC curves
	cv = StratifiedKFold(n_splits=5)
	classifier = LogisticRegression(solver='sag')

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

		print('roc_auc:', roc_auc)

		i += 1
		
	plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)

	cv_mean_auc = round(np.mean(aucs), 3)
	print('cv_mean_auc:', cv_mean_auc)


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

		
	pdf_filename = 'W' + str(win_len) + 'LR_ROC_gwRVIS.intergenic_vs_' + intol_class + discard_positive_str + '.pdf'
	fig.savefig(pdf_filename, bbox_inches='tight')
	
	return cv_mean_auc, intol_class_elems
	
	
	

	
	
	
	

if __name__ == '__main__':


	win_len = sys.argv[1]
	discard_positive = bool(int(sys.argv[2]))	# 1 or 0: discard or not windows with positive selection (gwRVIS > 0)
	

	intolerant_classes = ['ucne'] #, 'vista', 'utr', 'ccds', 'intron', 'lincrna']  
	print(intolerant_classes)

	
	
	# =============================== MAIN ANALYSIS ===============================	
	mean_auc_per_intol_class = {}
	
	for intol_class in intolerant_classes:

		input_classes = [intol_class, 'intergenic']
		obj = GwrvisDistributionPerClass(win_len=win_len, input_classes=input_classes)
		print('Sorted:', obj.sorted_genomic_classes)


		obj.read_gwrvis_scores_per_class(discard_positive_gwrvis=discard_positive)
		obj.compile_gwrvis_df()



		print('\nIntolerant class:', intol_class)
		toler_class = 'intergenic'

		
		
		cv_mean_auc, intol_class_elems = run_logit_regression_with_cv(obj.gwrvis_df, intol_class)
		
		print("Intergenic vs " + intol_class + ": " + str(cv_mean_auc))
		mean_auc_per_intol_class[intol_class] = [cv_mean_auc, intol_class_elems]

	
	
	discard_positive_str = ''
	if discard_positive:
		discard_positive_str = '.discard_positive'
	mean_aucs_out_file = 'mean_auc_per_intol_class.W' + str(win_len) + discard_positive_str + '.txt'
	
	
	with open(mean_aucs_out_file, 'w') as out_fh:
		out_fh.write('genomic_class\tauc\tclass_size\n')
		
		for k in sorted(mean_auc_per_intol_class.keys()):
			out_fh.write(k + '\t' + str(mean_auc_per_intol_class[k][0]) + '\t' + str(mean_auc_per_intol_class[k][1]) + '\n')
