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
from classifiers import Classifier

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
from custom_utils import create_out_dir



def prepare_feat_table_df(df, genomic_classes=[]):
	
	print('> All genomic classes:', df.genomic_class.unique())
	print('full_df:', df.shape)
	
	# filter for genomic classes if applicable
	if len(genomic_classes) > 0:
		df = df.loc[ df.genomic_class.isin(genomic_classes), :].copy()
	
	print('> Filtered genomic classes:', df.genomic_class.unique())
	print('subset_df:', df.shape)
	print('-----------')
	
	
	print(df.head())
	print(df.tail())
	print(df.columns)
	
	#df['uniq_id'] = df['chr'].astype(str) + '_' + df['start'].astype(str) + '_' + df['end'].astype(str)
	#print(len(df['uniq_id']))
	#print(len(df['uniq_id'].unique()))
	#df.to_csv('dbg_df.tsv', sep='\t', index=None)

	
	df[Y_label] = df[Y_label].astype(str).str.replace('Pathogenic.*', '1', regex=True)
	df[Y_label] = df[Y_label].astype(str).str.replace('Benign.*', '0', regex=True)
	df[Y_label] = df[Y_label].apply(pd.to_numeric, errors='coerce')

	return df
	
	
	
def prepare_extreme_intolerance_sets(df):

	df.sort_values(by=['gwrvis'], ascending=True, inplace=True)
	print(df.head())
	print(df.tail())
	

	
	



if __name__ == '__main__':

	config_file = sys.argv[1]	
	
	base_score = 'cadd' #'gwrvis'
	
	regression = False
	Y_label = 'clinvar_annot'
	if regression:
		Y_label = 'gwrvis'
	
	
	model_type = 'RandomForest' # 'RandomForest' # 'Logistic' #'RandomForest'
	include_vcf_extracted_features = False # default: False (including it for UTRs doesn't improve)
	use_only_base_score = False
	
	
	# If using a score other than gwRVIS, predict only using the score values without any other features 
	if base_score != 'gwrvis':
		use_only_base_score = True
	
	
	
	
	
	# =======  Dir Initialisation  =======
	out_dir = create_out_dir(config_file)    
	
	ml_data_dir = out_dir + '/ml_data'
	clinvar_feature_table_dir = ml_data_dir + '/clinvar_feature_tables'
	clinvar_ml_out_dir = ml_data_dir + '/clinvar-out'
	if not os.path.exists(clinvar_ml_out_dir):
		os.makedirs(clinvar_ml_out_dir)
	
	if base_score == 'gwrvis':
		full_feature_table_df = pd.read_csv(clinvar_feature_table_dir + '/full_feature_table.clinvar.bed', sep='\t', low_memory=False)
	else:
		full_feature_table_df = pd.read_csv(clinvar_feature_table_dir + '/full_feature_table.clinvar.' + base_score + '.bed', sep='\t', low_memory=False)

	print(full_feature_table_df.head())
	
		
	genomic_classes = ['utr', 'intergenic', 'lincrna', 'vista', 'ucne']
	#genomic_classes = ['intergenic']
	subset_df = prepare_feat_table_df(full_feature_table_df, genomic_classes=genomic_classes)

	
	
	# TODO:
	# - JARVIS performance per non-coding class vs any other score (w/ RandomForest)
	# - JARVIS performance across all coding classes compared to the other scores (still performs better :) )
	# - Add raw genomic sequence as a feature to JARVIS
	
	
	# - Subset genomic region within the Classifier class
	# - Create another class here to be called for different scores, receiving the results and then providing aggregate plots.
	#
	# Then:
	# - Train JARVIS with all HGMD pathogenic vs the ClinVar benign (or other set of benign variants). Then predict for all 3kb windows (with all the features already annotated) to rank them based on their probability score to be pathogenic.
	
	# Add HGMD as dataset
	# Add additional sets for control variants
	
	# Another project (almost):
	# - Predict most-intolerant vs most-tolerant from raw sequence only with CNNs (either as binary classification or regression).
	# The regression version may allow us to predict the gwRVIS score for regions that do not have variant data within a VCF.
	
	
	# ---------------------------------------------------------------------------------------	
	
	
	
	
	classifier = Classifier(Y_label, clinvar_ml_out_dir, regression=regression, model_type=model_type,
							include_vcf_extracted_features=include_vcf_extracted_features, 
							base_score=base_score, use_only_base_score=use_only_base_score)
	
	classifier.preprocess_data(subset_df)
	
	classifier.run_classification_with_cv()
	
	#classifier.train_and_predict()

