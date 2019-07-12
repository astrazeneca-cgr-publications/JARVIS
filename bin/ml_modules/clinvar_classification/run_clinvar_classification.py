import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_curve, auc, mean_squared_error
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
from custom_utils import create_out_dir


def prepare_feat_table_df(df, genomic_classes=[]):
	
	# filter for genomic classes if applicable
	if len(genomic_classes) > 0:
		print(genomic_classes)
		print(df.shape)
		df = df.loc[ df.genomic_class.isin(genomic_classes), :].copy()
		print(df.shape)
		print(df['genomic_class'].unique())
		print('-----------')
		
	cols_to_drop = ['chr', 'start', 'end', 'genomic_class']

	df.drop(cols_to_drop, axis=1, inplace=True)
	df[Y_label] = df[Y_label].astype(str).str.replace('Pathogenic.*', '1', regex=True)
	df[Y_label] = df[Y_label].astype(str).str.replace('Benign.*', '0', regex=True)
	df[Y_label] = df[Y_label].apply(pd.to_numeric, errors='coerce')

	return df
	
	
	
def prepare_extreme_intolerance_sets(df):

	df.sort_values(by=['gwrvis'], ascending=True, inplace=True)
	print(df.head())
	print(df.tail())
	


def run_logit_regression(df):

	feature_cols = [x for x in df.columns.values if x not in [Y_label, 'clinvar_annot'] ] #['gwrvis']
	print('Features:', feature_cols)
	X = df[feature_cols]
	y = df[Y_label].astype(int).values
	
	# logistic regression
	#model = LogisticRegression(C=1e9, solver='lbfgs', max_iter=10000)
	model = LinearRegression()
	
	model.fit(X, y)

	return model, X, y



def plot_roc_curve(model, df, X):

	print('RMSE:', mean_squared_error(df[Y_label], model.predict(X)))
	sys.exit()

	#df['pred_prob'] = model.predict_proba(X)[:, 1]

	fpr, tpr, thresholds = roc_curve(df[Y_label], df['pred_prob'])
	roc_auc = round(auc(fpr, tpr), 3)
	print("Area under the ROC curve : %f" % roc_auc)

	# Plot ROC curve
	fig, ax = plt.subplots(figsize=(10, 10))
	plt.plot(fpr, tpr, color='darkorange',
		 lw=2, label='ROC curve (area = %0.3f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic')
	plt.legend(loc="lower right", fontsize=14)
	plt.show()


	pdf_filename = clinvar_ml_out_dir + '/Logistic_Regression_ROC_gwRVIS.AUC_' + str(roc_auc) + '.pdf'

	fig.savefig(pdf_filename, bbox_inches='tight')




if __name__ == '__main__':

	config_file = sys.argv[1]	#"../config.yaml"
	
	Y_label = 'clinvar_annot'
	Y_label = 'gwrvis'
	
	
	out_dir = create_out_dir(config_file)    
	out_dir = '../../' + out_dir
	ml_data_dir = out_dir + '/ml_data'
	clinvar_feature_table_dir = ml_data_dir + '/clinvar_feature_tables'
	clinvar_ml_out_dir = ml_data_dir + '/clinvar-out'
	if not os.path.exists(clinvar_ml_out_dir):
		os.makedirs(clinvar_ml_out_dir)
	
	full_feature_table_df = pd.read_csv(clinvar_feature_table_dir + '/full_feature_table.clinvar.bed', sep='\t', low_memory=False)
		
	
	print(full_feature_table_df.head())
	print(full_feature_table_df.tail())
	
	print(full_feature_table_df['clinvar_annot'].unique())
	print('----------------------------\n')
	
	
	subset_df = prepare_feat_table_df(full_feature_table_df, genomic_classes=['utr', 'intergenic', 'lincrna', 'vista', 'ucne'])

	# TODO:
	# - Add standardisation
	# Add L1, L2 regulatisation
	#	- Look into feature selection
	# Run with random forest from sklearn
	#	- Look into feature importance with the default rf method
	
	# Explore linear regression to predict gwRVIS from the features (try avoiding circularity)
	
	model, X, y = run_logit_regression(subset_df)
	plot_roc_curve(model, subset_df, X)