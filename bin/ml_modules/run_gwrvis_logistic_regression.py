import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from custom_utils import create_out_dir



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


def run_logit_regression(df):

	feature_cols = ['gwrvis']
	X = df[feature_cols]
	y = df['gene_class']

	# logistic regression
	model = LogisticRegression(C=1e9, solver='lbfgs')
	model.fit(X, y)

	return model, X, y



def plot_roc_curve(model, df, X):

	df['pred_gene_class_prob'] = model.predict_proba(X)[:, 1]

	fpr, tpr, thresholds = roc_curve(df['gene_class'], df['pred_gene_class_prob'])
	roc_auc = round(auc(fpr, tpr), 2)
	print("Area under the ROC curve : %f" % roc_auc)

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
	discard_positive = int(sys.argv[2])	# 0 or 1
	
	out_dir = create_out_dir(config_file)    
	#full_genome_out_dir = '../' + out_dir + '/full_genome_out'
	full_genome_out_dir = out_dir + '/full_genome_out'


	full_gwrvis_df = pd.read_csv(full_genome_out_dir + '/Whole_genome_gwRVIS_per_class.csv')
	print(full_gwrvis_df.head())


	#print('Full Df size:', df.shape)
	#df = df.loc[ ~is_outlier(df.loc[:,'gwrvis']), :]
	#print('Df size (no outlier):', df.shape)


	intol_class = 'ucne'
	toler_class = 'intergenic'

	df = compile_gwrvis_differential_df(intol_class, toler_class)


	# ### Logistic Regression
	#df = df[ ~np.isnan(df.gwrvis) ]
	model, X, y = run_logit_regression(df)

	plot_roc_curve(model, df, X)
