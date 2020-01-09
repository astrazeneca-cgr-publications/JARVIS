import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import pickle

from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import roc_curve, auc, mean_squared_error
import sys
import os
os.environ['KMP_WARNINGS'] = 'off'
from classifiers import Classifier

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
from custom_utils import create_out_dir, get_config_params



class ClassificationWrapper:

	def __init__(self, config_file, base_score='gwrvis', model_type='DNN', genomic_classes=None, 
				Y_label='clinvar_annot', include_vcf_extracted_features=False, 
				exclude_base_score=False, filter_ccds_overlapping_variants=True):
		
		self.config_file = config_file
		self.base_score = base_score
		self.model_type = model_type
		self.genomic_classes = genomic_classes
		
		
		self.Y_label = Y_label		
		self.include_vcf_extracted_features = include_vcf_extracted_features
		self.exclude_base_score = exclude_base_score
		self.filter_ccds_overlapping_variants = filter_ccds_overlapping_variants
		
		self.use_only_base_score = True
		self.harmonise_options()
		self.init_dirs()
		
	

	def harmonise_options(self):
		
		if self.genomic_classes is None:
			self.genomic_classes = []
		
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

		
		self.base_out_models_dir = self.ml_data_dir + '/models'
		if not os.path.exists(self.base_out_models_dir): 			
			os.makedirs(self.base_out_models_dir)
	
	
	
	def read_input_data(self):
		
		if self.base_score in ['gwrvis', 'jarvis']:
			self.full_feature_table = pd.read_csv(self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.bed', sep='\t', low_memory=False)
		else:
			self.full_feature_table = pd.read_csv(self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.' + self.base_score + '.bed', sep='\t', low_memory=False)

		#print('>All features (prior to pre-processing):\n', self.full_feature_table.columns)
		#print(self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.bed')
		#print(self.full_feature_table.head())
	


	def subset_feat_table_df(self):
		
		print('> All genomic classes:', self.full_feature_table.genomic_class.unique())

		
		# Filter for genomic classes, if applicable
		if len(self.genomic_classes) > 0:
			self.df = self.full_feature_table.loc[ self.full_feature_table.genomic_class.isin(self.genomic_classes), :].copy()
		else:
			self.df = self.full_feature_table.copy()
		
		print('> Filtered genomic classes:', self.df.genomic_class.unique())
		
		# Drop NAs
		self.df.dropna(inplace=True)
		
		
		# Correct data types and convert Y-label strings to 1/0 values
		self.df[self.Y_label] = self.df[self.Y_label].astype(str).str.replace('Benign.*', '0', regex=True)
		self.df[self.Y_label] = self.df[self.Y_label].astype(str).str.replace('Pathogenic.*', '1', regex=True)
			
	
		if self.base_score in ['gwrvis', 'jarvis']:
			cur_full_feature_file = self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.' + '_'.join(self.genomic_classes) + '.bed'
			clean_feature_file = self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.' + '_'.join(self.genomic_classes) + '.clean.bed'

			# Get the blacklisted regions (gwRVIS windows which contain both CCDS and non-coding variants)
			os.system("./jarvis/variant_classification/get_overlapping_variant_windows.sh " + genomic_classes_log + " " + self.out_dir + " " + self.clinvar_feature_table_dir + " " + patho_benign_sets)

		else:
			cur_full_feature_file = self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.' + self.base_score + "." + '_'.join(self.	genomic_classes) + '.bed'
			clean_feature_file = self.clinvar_feature_table_dir + '/full_feature_table.' + patho_benign_sets + '.' + self.base_score + "." + '_'.join(self.genomic_classes) + '.clean.bed'
			
		self.df.to_csv(cur_full_feature_file, sep='\t', index=False, header=False)

		
		
		if self.filter_ccds_overlapping_variants:
			# Filter-out the blacklisted regions (gwRVIS windows which contain both CCDS and non-coding variants)
			clean_file_is_initialised = False
			
			if 'intergenic' in self.genomic_classes:

				os.system("bedtools subtract -a " + cur_full_feature_file + " -b " + self.clinvar_feature_table_dir + "/blacklist.intergenic_overlaping_ccds.bed > " + clean_feature_file)
				clean_file_is_initialised = True
				
			if 'utr' in self.genomic_classes:
				if clean_file_is_initialised:
					os.system("bedtools subtract -a " + clean_feature_file + " -b " + self.clinvar_feature_table_dir + "/blacklist.utr_overlaping_ccds.bed > " + clean_feature_file + '.tmp')
					
					os.system("mv " + clean_feature_file + ".tmp " + clean_feature_file)
				else:
					os.system("bedtools subtract -a " + cur_full_feature_file + " -b " + self.clinvar_feature_table_dir + "/blacklist.utr_overlaping_ccds.bed > " + clean_feature_file)
					
		
		
			# Replace df with clean version only when intergenic or UTR are in the genomic_classes
			if 'intergenic' in self.genomic_classes or 'utr' in self.genomic_classes:
				original_columns = self.df.columns
			
				self.df = pd.read_csv(clean_feature_file, sep='\t', header=None)
				self.df.columns = original_columns
		
		print(self.df.head())
		print(self.df.tail())
		print(self.df.shape)
		

		
	
	def run_classifier(self):
	
		classifier_out_dir = self.clinvar_ml_out_dir + '/' + '_'.join(self.genomic_classes)
		if not os.path.exists(classifier_out_dir):
			os.makedirs(classifier_out_dir)

		# @anchor-1
		#out_models_dir = self.base_out_models_dir + '/' + 'intergenic_utr_lincrna_ucne_vista'  #+ '_'.join(self.genomic_classes)
		out_models_dir = self.base_out_models_dir + '/' + 'intergenic'  #+ '_'.join(self.genomic_classes)
		"""
		    -- JARVIS:
		       > "intergenic_utr_lincrna_ucne_vista" is the best model for 'utr'--D3000-struct (0.675)
		       > "intergenic_utr_lincrna_ucne_vista" is the best model for 'intergenic,utr'--D3000-struct (0.649)
		       > "ccds" is (probably) the best model for 'ccds'--D3000-struct (0.565)

		"""
		if not os.path.exists(out_models_dir): 			
			os.makedirs(out_models_dir)
			
	
		classifier = Classifier(self.Y_label, classifier_out_dir, out_models_dir, base_score=self.base_score,
							model_type=self.model_type,
							use_only_base_score=self.use_only_base_score,
							include_vcf_extracted_features=self.include_vcf_extracted_features, 
							exclude_base_score=self.exclude_base_score)
		
		classifier.preprocess_data(self.df)
		
		classifier.init_model()
		
		
		classifier.run_classification_with_cv()
		
		self.score_print_name = classifier.score_print_name
		self.mean_tpr = classifier.mean_tpr
		self.mean_fpr = classifier.mean_fpr	
		self.mean_auc = classifier.mean_auc
		
		self.metrics_list = classifier.metrics_list
		

	
		
	def run(self):
	
		if self.Y_label == 'clinvar_annot':
			# Clinvar-based classification
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
	

	pdf_filename = clinvar_ml_out_dir + '/' + '_'.join(genomic_classes) + '.all_scores_classification.' + Y_label + '.pdf'
	fig.savefig(pdf_filename, bbox_inches='tight')
	

	
	
def check_and_save_performance_metrics(metrics_per_score, genomic_classes, clinvar_feature_table_dir):

	for score, metrics in metrics_per_score.items():
		print('\n>', score)
		#print('Metrics:', metrics)
	
		#print("metrics[0].keys():", list(metrics[0].keys()))
		
		
		for metric in metrics[0].keys():
			#print('Metric:', metric)
			cur_metric_list = []
			
			for sublist in metrics:
				cur_metric_list.append(sublist[metric])
				
				cur_metric_list = [x if x != 'NA' else 0.5 for x in cur_metric_list]
				
				
			#print(cur_metric_list)
			print('Avg.', metric, ':', np.mean(cur_metric_list))
				
	metrics_out_file = clinvar_feature_table_dir + '/all_but_dnn_performance_metrics.' + '_'.join(genomic_classes) + '.pkl'
	with open(metrics_out_file, 'wb') as handle:
		pickle.dump(metrics_per_score, handle, protocol=pickle.HIGHEST_PROTOCOL)
		
	print('Metrics saved at:', metrics_out_file)
		
	
	
	

if __name__ == '__main__':

	config_file = sys.argv[1]
	filter_ccds_overlapping_variants = bool(int(sys.argv[2]))
	model_type = sys.argv[3] #'DNN' (default) # 'RF (RandomForest)', 'Logistic', 'DNN'

	run_params = get_config_params(config_file)
	
	genomic_classes_log = run_params['genomic_classes']
	Y_label = run_params['Y_label']
	
	pathogenic_set = run_params['pathogenic_set']
	benign_set = run_params['benign_set']
	print('Pathogenic set: ' + pathogenic_set)
	print('Benign set: ' + benign_set)
	patho_benign_sets = pathogenic_set + '_' + benign_set


	if Y_label == 'clinvar_annot':
		#genomic_classes_lists =  [ ['intergenic'], ['utr'], ['lincrna'], ['intergenic', 'utr', 'lincrna', 'ucne', 'vista'], ['intergenic', 'utr', 'lincrna', 'ucne', 'vista', 'intron'], ['intergenic', 'utr', 'lincrna', 'ucne', 'vista', 'ccds', 'intron'], ['ccds'], ['intron'] ] 
		# @anchor-2
		genomic_classes_lists =  [ ['intergenic'] ]  # TEMP

	
	
	hg_version = run_params['hg_version']
	if hg_version == 'hg19':
		#all_base_scores = ['ncER_10bp', 'cdts', 'linsight', 'gwrvis', 'jarvis', 'cadd', 'dann', 'phyloP46way', 'phastCons46way', 'orion'] 
		# @anchor-3
		all_base_scores = ['jarvis'] #['jarvis']  # TEMP
	else:
		all_base_scores = ['gwrvis', 'jarvis']

	# Ad-hoc: exclude 'lincrna' when running for ncRVIS
	#genomic_classes_lists =  [ ['intergenic'], ['utr'], ['utr', 'intergenic', 'lincrna', 'vista', 'ucne'], ['ccds'], ['intron'] ] 
	#all_base_scores = ['ncRVIS'] 
	

	# include_vcf_extracted_features -- default: False (including it for UTRs doesn't improve)

	
	
	
	
	
	for genomic_classes in genomic_classes_lists:
		print('\n**********************************************************************\n>>\t\t\t\t' + ' '.join(genomic_classes) + '\n**********************************************************************\n')
		score_list, fpr_list, tpr_list, auc_list = [], [], [], []
		
		
		metrics_per_score = {}
		
		for base_score in all_base_scores:

			print('>>>>>>>  ' + base_score + '\n')

			#try:  # 16 lines in try
			clf_wrapper = ClassificationWrapper(config_file, base_score=base_score, model_type=model_type, 
												genomic_classes=genomic_classes,
												Y_label=Y_label, 
												include_vcf_extracted_features=False, 
												exclude_base_score=False,
												filter_ccds_overlapping_variants=filter_ccds_overlapping_variants)
												
			clf_wrapper.run()
		
			if clf_wrapper.mean_auc is not None:
				score_list.append(clf_wrapper.score_print_name)
				fpr_list.append(clf_wrapper.mean_fpr)
				tpr_list.append(clf_wrapper.mean_tpr)
				auc_list.append(clf_wrapper.mean_auc)

			metrics_per_score[base_score] = clf_wrapper.metrics_list
			#except:
			#	print("\n\n[Exception] in " + ','.join(genomic_classes) + " for score: " + base_score + "\n") 
		
		plot_roc_curve(score_list, fpr_list, tpr_list, auc_list, genomic_classes, clf_wrapper.clinvar_ml_out_dir, all_base_scores)
		
		check_and_save_performance_metrics(metrics_per_score, genomic_classes, clf_wrapper.clinvar_feature_table_dir)

		
