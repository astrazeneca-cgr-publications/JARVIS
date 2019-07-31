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

	config_file = sys.argv[1]	#"../../config.yaml"
	
	base_score = 'gwrvis' #'gwrvis'
	
	regression = False
	Y_label = 'clinvar_annot'
	if regression:
		Y_label = 'gwrvis'
	
	
	out_dir = create_out_dir(config_file)    
	out_dir = out_dir
	#out_dir = '../../' + out_dir
	
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
	
		
	#genomic_classes = ['utr', 'intergenic', 'lincrna', 'vista', 'ucne']
	genomic_classes = ['intergenic']
	subset_df = prepare_feat_table_df(full_feature_table_df, genomic_classes=genomic_classes)

	
	
	# ---------------------------------------------------------------------------------------	
	
	
	model_type = 'RandomForest' # 'Logistic' #'RandomForest'
	include_vcf_extracted_features = False
	use_only_base_score = False
	
	
	classifier = Classifier(Y_label, clinvar_ml_out_dir, regression=regression, model_type=model_type,
							include_vcf_extracted_features=include_vcf_extracted_features, 
							base_score=base_score, use_only_base_score=use_only_base_score)
	classifier.split_into_train_test_sets(subset_df)
	
	classifier.train_and_predict()
