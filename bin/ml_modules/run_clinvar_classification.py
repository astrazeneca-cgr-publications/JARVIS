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


def prepare_extreme_intolerance_sets(df):

	df.sort_values(by=['gwrvis'], ascending=True, inplace=True)
	print(df.head())
	print(df.tail())


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
	
	out_dir = create_out_dir(config_file)    
	out_dir = '../' + out_dir
	ml_data_dir = out_dir + '/ml_data'
	feature_table_dir = ml_data_dir + '/feature_tables'
	
	full_gwrvis_feature_table_df = pd.read_csv(feature_table_dir + '/full_gwrvis_and_regulatory_features.All_genomic_classes.tsv', sep='\t', low_memory=False)
	#full_gwrvis_feature_table_df = full_gwrvis_feature_table_df[np.isfinite(full_gwrvis_feature_table_df['gwrvis'])]
	#full_gwrvis_feature_table_df['gwrvis'] = pd.to_numeric(full_gwrvis_feature_table_df['gwrvis'])
	
	print(full_gwrvis_feature_table_df.shape)
	full_gwrvis_feature_table_df.dropna(inplace=True)
	print(full_gwrvis_feature_table_df.shape)
	
	print(full_gwrvis_feature_table_df.head())
	print(full_gwrvis_feature_table_df.tail())
	print(full_gwrvis_feature_table_df.info())
	print('----------------------------\n')
	
	prepare_extreme_intolerance_sets(full_gwrvis_feature_table_df)
	sys.exit()


	intol_class = 'ucne'
	toler_class = 'intergenic'
	
	df = compile_gwrvis_differential_df(intol_class, toler_class)


	# ### Logistic Regression
	#df = df[ ~np.isnan(df.gwrvis) ]
	model, X, y = run_logit_regression(df)

	plot_roc_curve(model, df, X)
