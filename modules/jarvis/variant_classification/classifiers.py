import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from collections import Counter
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import load_model

from scipy import interp, stats
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV, LinearRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from imblearn.under_sampling import RandomUnderSampler

import warnings
warnings.filterwarnings("error")
import pickle
import sys
import os



sys.path.insert(1, os.path.join(sys.path[0], '../..'))
from custom_utils import create_out_dir
	

	

class Classifier:
	
	def __init__(self, Y_label, out_dir, out_models_dir, base_score='gwrvis', model_type='RF', 
				 use_only_base_score=True, include_vcf_extracted_features=False, exclude_base_score=False, 
				 use_pathogenicity_trained_model=False, use_conservation_trained_model=False,
				 predict_on_test_set=False):
				 
		self.out_dir = out_dir
		self.out_models_dir = out_models_dir
		self.Y_label = Y_label
		self.model_type = model_type
		self.include_vcf_extracted_features = include_vcf_extracted_features
		self.base_score = base_score
		self.use_only_base_score = use_only_base_score
		self.exclude_base_score = exclude_base_score
		self.use_pathogenicity_trained_model = use_pathogenicity_trained_model
		self.use_conservation_trained_model = use_conservation_trained_model
		self.predict_on_test_set = predict_on_test_set
		

			
	
	def preprocess_data(self, df):
	
	
		if self.use_only_base_score:		
			df.dropna(inplace=True)
			self.feature_cols = [self.base_score]
		
		else:
			#cols_to_drop = ['chr', 'start', 'end', 'genomic_class', self.Y_label, 'clinvar_annot']

			cols_to_drop = ['chr', 'start', 'end', 'genomic_class', self.Y_label, 'clinvar_annot']
			

	
			vcf_dependent_cols = ['common_variants', 'common_vs_all_variants_ratio', 
								  'all_variants', 'mean_ac', 'mean_af', 'bin_1', 
								  'bin_2', 'bin_3', 'bin_4', 'bin_5', 'bin_6']
								  
			if not self.include_vcf_extracted_features:
				cols_to_drop.extend(vcf_dependent_cols)
				
			if self.exclude_base_score:
				cols_to_drop.append('gwrvis')
				
			self.feature_cols = [x for x in df.columns.values if x not in cols_to_drop ]
		


		#print('Features:', self.feature_cols)
		#print('Label:', self.Y_label)
	
		print(Counter(df[self.Y_label]))
	
		# BETA
		# Retaining only most important features (improves AUC only by +0.001)
		#self.feature_cols = ['gwrvis', 'gc_content', 'cpg', 'mut_rate', 'cpg_islands', 'H3K4me2']
			
		#print(df.loc[ df[self.Y_label] == '1', :])
		

		
		# ==== BETA - Temporary (Subset variants to those used by Linsight) ====
		"""
		df['aux'] = df['chr'] + '_' + df['start'].astype(str) + '_' + df['end'].astype(str)
		
		if self.base_score == 'jarvis':

			tt = pd.read_csv('linsight.all_non_coding.coords', index_col=0, header=None)
			tt.columns = ['aux']
			print(df.shape)

			df = df.loc[ df.aux.isin(tt.aux.values), : ].copy()
			print(df.head())
			print(df.shape)
			print(tt.head())
			print(tt.shape)

	
		if self.base_score == 'linsight':
			df['aux'].to_csv('linsight.all_non_coding.coords', header=False)
		
		df.drop(['aux'], axis=1, inplace=True)
		"""
		# =======================================================================



		self.X = df[self.feature_cols].values
		self.y = df[self.Y_label].astype(int).values

		print(self.X[0])
		# Standardise features -- [Deprecated - distorts feature importance and decreases AUC performance]
		#self.X = stats.zscore(self.X, axis=1)
		
		# Normalise features
		min_max_scaler = preprocessing.MinMaxScaler()
		self.X = min_max_scaler.fit_transform(self.X)

		print(self.X[0])

		
		
		# Fix class imbalance (with over/under-sampling minority/majority class)
		positive_set_size = (self.y == 1).sum()
		negative_set_size = (self.y == 0).sum()
		pos_neg_ratio = 1/1

		
		
		if (positive_set_size / negative_set_size < pos_neg_ratio) or (negative_set_size / positive_set_size < pos_neg_ratio):
			print('\n> Fixing class imbalance ...')
			print('Imbalanced sets: ', sorted(Counter(self.y).items()))
			rus = RandomUnderSampler(random_state=0, sampling_strategy=pos_neg_ratio)
			self.X, self.y = rus.fit_resample(self.X, self.y)
			print('Balanced sets:', sorted(Counter(self.y).items()))
				
				
		
	
	def init_model(self):
		"""
			Initialise classifier model based on input base_score and model_type
		"""
			
		#rf_params = dict(n_estimators=200, max_depth=3, random_state=0)
		rf_params = dict(n_estimators=100, max_features=5, max_depth=2)
		gb_params = dict(n_estimators=100, max_features=5, max_depth=2)
	
	
		# Use logistic regression for all scores except for JARVIS
		if self.base_score != 'jarvis':
			self.model_type = 'Logistic'
			self.model = LogisticRegression(C=1, solver='lbfgs', max_iter=10000)
		else:
			if self.model_type == 'RF':
				self.model = RandomForestClassifier(**rf_params)	
				
			elif self.model_type == 'GB':
				self.model = GradientBoostingClassifier(**gb_params)	

			elif self.model_type == 'Logistic':
				self.model = LogisticRegression(C=1e9, solver='lbfgs', max_iter=10000)
			
		print('Model:', self.model_type, '\n')
		
		
		
		self.score_print_name = self.base_score
		
		if self.base_score == 'gwrvis':
			self.score_print_name = 'gwRVIS'
						
		if self.base_score == 'cadd':
			self.score_print_name = self.base_score.upper()
		
		if self.base_score == 'orion':
			self.score_print_name = self.base_score.capitalize()
			
		if self.base_score == 'jarvis':
			self.base_score = 'gwrvis'
			self.score_print_name = 'JARVIS'
	
	
	
	
		
		
	def run_classification_with_cv(self, cv_splits=5, cv_repeats=1):
		
		cv = StratifiedKFold(n_splits=cv_splits)
		
		
		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)
		metrics_list = []
		

		y_label_lists = []
		y_proba_lists = []
	

		fig, ax = plt.subplots(figsize=(10, 10))


		if self.base_score == 'gwrvis':
			print(self.X.shape)
			print('Features:', self.feature_cols)
			print('================================================')



	
	
		# n-repeated CV
		for n in range(cv_repeats):
			print('CV - Repeat:', str(n+1))
			fold = 0
			for train, test in cv.split(self.X, self.y):

			
				if self.base_score not in ['gwrvis', 'jarvis']:
					probas_ = self.model.fit(self.X[train], self.y[train]).predict_proba(self.X[test])
					
				else:
					if self.use_pathogenicity_trained_model:
							
							# -- At the moment use RF-based model for JARVIS (and LR for gwRVIS anyway)
							#model_out_file = self.out_models_dir + '/' + self.score_print_name + '-' + self.model_type + '.model'
							full_out_models_dir = "../out/topmed-NEW_ClinVar_pathogenic-denovodb_benign-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/models/intergenic_utr_lincrna_ucne_vista"
							
							model_out_file = full_out_models_dir + '/' + self.score_print_name + '-' + self.model_type + '.model'
							print("\n>> Loading PATHOGENICITY-trained model from file:", model_out_file)
							
							# scikit-learn model
							with open(model_out_file, 'rb') as fh:     
								self.model = pickle.load(fh)	
							probas_ = self.model.predict_proba(self.X[test])
						
					else:
					
					
						if self.predict_on_test_set:
						
							# Predict on full X dataset, without cross-validation splits
							print("Predicting on test set...")
							#full_out_models_dir = self.out_dir + '/../../models/intergenic_utr_lincrna_ucne_vista/' 
							
							#full_out_models_dir = "/projects/cgr/users/kclc950/JARVIS/out/topmed-ClinVar_pathogenic-denovodb_benign-winlen_1000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/models/intron_intergenic_utr_lincrna_ucne_vista/"
							full_out_models_dir = "../out/topmed-NEW_ClinVar_pathogenic-denovodb_benign-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/models/intergenic_utr_lincrna_ucne_vista"
							
							model_out_file = full_out_models_dir + '/' + self.score_print_name + '-' + self.model_type + '.model'
							print('model_out_file:', model_out_file)

							with open(model_out_file, 'rb') as fh:
								self.model = pickle.load(fh)

							probas_ = self.model.predict_proba(self.X)
							
							
							# Compute ROC curve and area the curve
							fpr, tpr, thresholds = roc_curve(self.y, probas_[:, 1])
								
							
							tprs.append(interp(mean_fpr, fpr, tpr))
							tprs[-1][0] = 0.0
							roc_auc = round(auc(fpr, tpr), 3)
							aucs.append(roc_auc)
							plt.plot(fpr, tpr, lw=1, alpha=0.3,
									 label='ROC fold %d (AUC = %0.2f)' % (fold, roc_auc))

							
							# Evaluate predictions on test and get performance metrics
							metrics = self.test_and_evaluate_model(probas_, self.y)
							metrics['auc'] = roc_auc
							metrics_list.append(metrics)		 
						 
							break
					
						else:
							probas_ = self.model.fit(self.X[train], self.y[train]).predict_proba(self.X[test])
						



				y_label_lists.append(self.y[test].tolist())
				y_proba_lists.append(probas_[:, 1].tolist())
	

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
				
		
				if self.predict_on_test_set:
					break
		

		self.y_label_lists = y_label_lists
		self.y_proba_lists = y_proba_lists


		
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
		plt.title(self.score_print_name + ': ' + str(cv_splits) + '-fold Cross-Validation ROC Curve')
		plt.legend(loc="lower right")
		plt.show()
		
		
		pdf_filename = self.out_dir + '/' + self.model_type + '.' + self.Y_label
						
		if self.include_vcf_extracted_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')
		plt.close('all')
		
		
		print('Mean AUC:', self.mean_auc)
		
		
		if self.model_type in ['RF', 'GB', 'Logistic']:
			self.get_feature_importances()

		self.mean_tpr = mean_tpr
		self.mean_fpr = mean_fpr
		print('=====================\n\n')
	
		self.metrics_list = metrics_list
		#print("\nMetrics: ", metrics_list)
	
	
		# =========== TRAIN FULL MODEL FOR JARVIS and gwRVIS ===========
		if self.base_score in ['gwrvis', 'jarvis'] and not (self.use_conservation_trained_model or self.use_pathogenicity_trained_model or self.predict_on_test_set):
		
			# Save the model only when a pre-trained one has NOT been defined to be used for prediction
			self.model.fit(self.X, self.y)
			model_out_file = self.out_models_dir + '/' + self.score_print_name + '-' + self.model_type + '.model'
			with open(model_out_file, 'wb') as fh:
				pickle.dump(self.model, fh)

			print("Saved full model into:", model_out_file, '\n')
			
			
	
	
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
	
	
		
	
	def get_feature_importances(self):
	
		print("\n> Feature importances:")
		if self.model_type == 'Logistic':
			importances = self.model.coef_.reshape(-1, 1)[:, 0]
		elif self.model_type in ['RF', 'GB']:
			importances = self.model.feature_importances_
		
		
		importances_series = pd.Series(importances, index=self.feature_cols)
		importances_series.sort_values(inplace=True)
		print('Feature importance:\n', importances_series)
		
		fig, ax = plt.subplots(figsize=(14, 10))
		
		importances_series.sort_values(ascending=True, inplace=True) 

		importances_series.plot.barh()
		plt.show()

		pdf_filename = self.out_dir + '/' + self.model_type + '_importance_scores.' + self.score_print_name + '.' + self.model_type + '.' + self.Y_label

		if self.include_vcf_extracted_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')

