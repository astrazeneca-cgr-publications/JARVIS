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



class ClassificationWrapper:

	def __init__(self, config_file, base_score='gwrvis', model_type='RandomForest', genomic_classes=None, 
				Y_label='clinvar_annot', use_only_base_score=True, include_vcf_extracted_features=False, regression=False):
		
		self.config_file = config_file
		self.base_score = base_score
		self.model_type = model_type
		self.genomic_classes = genomic_classes
		
		
		self.Y_label = Y_label		
		self.use_only_base_score = use_only_base_score
		self.include_vcf_extracted_features = include_vcf_extracted_features
		self.regression = regression

		
		self.harmonise_options()
		self.init_dirs()
		
	

	def	harmonise_options(self):
		
		if self.genomic_classes is None:
			self.genomic_classes = []
		
		# If regression == True, then the target variable is 'gwrvis'
		if self.regression:
			self.Y_label = 'gwrvis'
		else:
			self.Y_label = 'clinvar_annot' # Default: for classification of pathogenic-vs-benign ClinVar variants
	
		# If using a score other than gwRVIS (e.g. CADD, phylop, etc.) 
		# predict using only the score's values without any other features 
		if self.base_score == 'jarvis':
			self.use_only_base_score = False
	
	
	
	def init_dirs(self):
		""" 
			Dir Initialisation
		"""
		
		self.out_dir = create_out_dir(self.config_file, create_dirs=False)    
		self.ml_data_dir = self.out_dir + '/ml_data'
		
		self.clinvar_feature_table_dir = self.ml_data_dir + '/clinvar_feature_tables'
		
		self.clinvar_ml_out_dir = self.ml_data_dir + '/clinvar-out'
		if not os.path.exists(self.clinvar_ml_out_dir):
			os.makedirs(self.clinvar_ml_out_dir)
	
	
	
	def read_input_data(self):
		
		if self.base_score in ['gwrvis', 'jarvis']:
			self.full_feature_table = pd.read_csv(self.clinvar_feature_table_dir + '/full_feature_table.clinvar.bed', sep='\t', low_memory=False)
		else:
			self.full_feature_table = pd.read_csv(self.clinvar_feature_table_dir + '/full_feature_table.clinvar.' + base_score + '.bed', sep='\t', low_memory=False)

		print('>All features (prior to pre-processing):\n', self.full_feature_table.columns)
		
	

	def subset_feat_table_df(self):
		
		print('> All genomic classes:', self.full_feature_table.genomic_class.unique())

		
		# Filter for genomic classes, if applicable
		if len(self.genomic_classes) > 0:
			self.df = self.full_feature_table.loc[ self.full_feature_table.genomic_class.isin(self.genomic_classes), :].copy()
		else:
			self.df = self.full_feature_table.copy()
			
		print('> Filtered genomic classes:', self.df.genomic_class.unique())
		
		
		# Correct data types and convert Y-label strings to 1/0 values
		self.df[self.Y_label] = self.df[self.Y_label].astype(str).str.replace('Pathogenic.*', '1', regex=True)
		self.df[self.Y_label] = self.df[self.Y_label].astype(str).str.replace('Benign.*', '0', regex=True)
		self.df[self.Y_label] = self.df[self.Y_label].apply(pd.to_numeric, errors='coerce')
	
	
	
	def run_classifier(self):
	
		classifier_out_dir = self.clinvar_ml_out_dir + '/' + '_'.join(self.genomic_classes)
		if not os.path.exists(classifier_out_dir):
			os.makedirs(classifier_out_dir)
			
	
		classifier = Classifier(self.Y_label, classifier_out_dir, base_score=self.base_score,
							model_type=self.model_type,
							use_only_base_score=self.use_only_base_score,
							include_vcf_extracted_features=self.include_vcf_extracted_features, 
							regression=self.regression)
		
		classifier.preprocess_data(self.df)
		
		classifier.run_classification_with_cv()
		
		self.score_print_name = classifier.score_print_name
		self.mean_tpr = classifier.mean_tpr
		self.mean_fpr = classifier.mean_fpr	
		self.mean_auc = classifier.mean_auc

		
	def run(self):
	
		self.read_input_data()
		
		self.subset_feat_table_df()
		
		self.run_classifier()
		
		
		
		
		
		
def plot_roc_curve(score_list, fpr_list, tpr_list, auc_list, genomic_classes, clinvar_ml_out_dir, all_base_scores):

	rvis_colors = ['#e31a1c', '#4292c6']
	colors = [c for c in sns.color_palette("Paired", 12).as_hex() if c not in rvis_colors]	
	class_colors = dict( zip(score_list, colors[:len(score_list)]) )
	class_colors['gwRVIS'] = '#e31a1c'
	class_colors['JARVIS'] = '#4292c6'
		
		
	# Plot ROC curve
	fig, ax = plt.subplots(figsize=(10, 10))


	for i in range(len(score_list)):
		cur_score = score_list[i]
		fpr = fpr_list[i]
		tpr = tpr_list[i]
		roc_auc = auc_list[i]
			
		lw = 1
		if cur_score == 'gwRVIS':
			lw = 1.5
		elif cur_score == 'JARVIS':
			lw = 2
		
		plt.plot(fpr, tpr, color=class_colors[cur_score],
				 lw=lw, label=cur_score + ' (AUC = %0.3f)' % roc_auc)


	plt.plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('ROC Curves - ' + ';'.join(genomic_classes))
	plt.legend(loc="lower right", fontsize=14)
	plt.show()
	

	pdf_filename = clinvar_ml_out_dir + '/' + '_'.join(genomic_classes) + '.all_scores_classification.pdf'
	fig.savefig(pdf_filename, bbox_inches='tight')
	
	
	
	
	
	
	
	
	

if __name__ == '__main__':

	config_file = sys.argv[1]	
		
	genomic_classes_lists = [ ['intergenic'], ['utr'], ['utr', 'intergenic', 'lincrna', 'vista', 'ucne'] ]
	all_base_scores = ['jarvis', 'gwrvis', 'cadd', 'phyloP46way', 'phastCons46way', 'orion']
	
	# include_vcf_extracted_features -- default: False (including it for UTRs doesn't improve)
	# use_only_base_score -- default: False (relevant only for gwRVIS; 'use_only_base_score' is always True for all other scores
	# regression -- default: False
	

	model_type = 'RandomForest' # 'RandomForest' (default), 'Logistic' 
	
	
	for genomic_classes in genomic_classes_lists:
	
		score_list, fpr_list, tpr_list, auc_list = [], [], [], []
		
		for base_score in all_base_scores:
	
			clf_wrapper = ClassificationWrapper(config_file, base_score=base_score, model_type=model_type, 
												genomic_classes=genomic_classes,
												Y_label='clinvar_annot', 
												use_only_base_score=True, 
												include_vcf_extracted_features=False, 
												regression=False)
												
			clf_wrapper.run()
		
			score_list.append(clf_wrapper.score_print_name)
			fpr_list.append(clf_wrapper.mean_fpr)
			tpr_list.append(clf_wrapper.mean_tpr)
			auc_list.append(clf_wrapper.mean_auc)

			
		print(fpr_list)
		print(tpr_list)
		print(auc_list)
		
		
		
		plot_roc_curve(score_list, fpr_list, tpr_list, auc_list, genomic_classes, clf_wrapper.clinvar_ml_out_dir, all_base_scores)
		

		
		
	
	
	# TODO:
	# - Add raw genomic sequence as a feature to JARVIS
		# > Need to add DNN as classifier in that case
		
	
	# - JARVIS performance per non-coding class vs any other score (w/ RandomForest)  [DONE]
	# - JARVIS performance across all non-coding classes compared to the other scores (still performs better :) ) [DONE]
	# - Subset genomic region within the Classifier class [DONE]
	# - Create another class here to be called for different scores, receiving the results and then providing aggregate plots. [DONE]
	
	
	# Then:
	# - Train JARVIS with all HGMD pathogenic vs the ClinVar benign (or other set of benign variants). Then predict for all 3kb windows (with all the features already annotated) to rank them based on their probability score to be pathogenic.
	
	# Add HGMD as dataset
	# Add additional sets for control variants
	# Add annotation for Histon marks/Methylation from other cell types too.
	
	# Another project (almost):
	# - Predict most-intolerant vs most-tolerant from raw sequence only with CNNs (either as binary classification or regression).
	# The regression version may allow us to predict the gwRVIS score for regions that do not have variant data within a VCF.
	# ---------------------------------------------------------------------------------------	