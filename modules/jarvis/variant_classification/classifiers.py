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

from scipy import interp
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV, LinearRegression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler

import warnings
warnings.filterwarnings("error")
import pickle
import sys
import os

#import tensorflow as tf


sys.path.insert(1, os.path.join(sys.path[0], '../..'))
from custom_utils import create_out_dir
	

	

class Classifier:
	
	def __init__(self, Y_label, out_dir, out_models_dir, base_score='gwrvis', model_type='RF', 
				 use_only_base_score=True, include_vcf_extracted_features=False, exclude_base_score=False):
				 
		self.out_dir = out_dir
		self.out_models_dir = out_models_dir
		self.Y_label = Y_label
		self.model_type = model_type
		self.include_vcf_extracted_features = include_vcf_extracted_features
		self.base_score = base_score
		self.use_only_base_score = use_only_base_score
		self.exclude_base_score = exclude_base_score

			
	
	def preprocess_data(self, df):
	
		if self.use_only_base_score:		
			df.dropna(inplace=True)
			self.feature_cols = [self.base_score]
		
		else:
			cols_to_drop = ['chr', 'start', 'end', 'genomic_class', self.Y_label, 'clinvar_annot']
	
			vcf_dependent_cols = ['common_variants', 'common_vs_all_variants_ratio', 
								  'all_variants', 'mean_ac', 'mean_af', 'bin_1', 
								  'bin_2', 'bin_3', 'bin_4', 'bin_5', 'bin_6']
								  
			if not self.include_vcf_extracted_features:
				cols_to_drop.extend(vcf_dependent_cols)
				
			if self.exclude_base_score:
				cols_to_drop.append('gwrvis')
				
			self.feature_cols = [x for x in df.columns.values if x not in cols_to_drop ]
		
		print('Features:', self.feature_cols)
		print('Label:', self.Y_label)
	
		print(Counter(df[self.Y_label]))
	
		# BETA
		# Retaining only most important features (improves AUC only by +0.001)
		#self.feature_cols = ['gwrvis', 'gc_content', 'cpg', 'mut_rate', 'cpg_islands', 'H3K4me2']
			
		self.X = df[self.feature_cols].values
		self.y = df[self.Y_label].astype(int).values

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
			
		rf_params = dict(n_estimators=100, max_depth=2, random_state=0)
	
	
		# Use logistic regression for all scores except for JARVIS
		if self.base_score != 'jarvis':
			self.model_type = 'Logistic'
			self.model = LogisticRegression(C=1, solver='lbfgs', max_iter=10000)
		else:
			if self.model_type == 'RF':
				self.model = RandomForestClassifier(**rf_params)	
				
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
	
	
	
	
		
		
	def run_classification_with_cv(self, cv_splits=5, cv_repeats=5):
		
		cv = StratifiedKFold(n_splits=cv_splits)
		
		
		tprs = []
		aucs = []
		mean_fpr = np.linspace(0, 1, 100)
		metrics_list = []
		
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

			
				probas_ = self.model.fit(self.X[train], self.y[train]).predict_proba(self.X[test])

				# BETA
				"""
				if self.base_score == 'gwrvis':
					# TODO:
					self.file_annot = 'D3000.no_zeros'

					model_out_file = self.out_models_dir + '/' + self.score_print_name + '-' + self.model_type + '.' + self.file_annot + '.model'
					
					#print("-- Loading pre-built model from file:", model_out_file)
					with open(model_out_file, 'rb') as fh:     
						self.model = pickle.load(fh)	
					probas_ = self.model.predict_proba(self.X[test])
				"""

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
		plt.close()
		
		
		pdf_filename = self.out_dir + '/' + self.model_type + '_ROC.' + self.score_print_name + \
						'.AUC_' + str(self.mean_auc) + '.' + self.Y_label
						
		if self.include_vcf_extracted_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')
		
		
		print('Mean AUC:', self.mean_auc)
		
		
		if self.model_type == 'RF':
			self.get_feature_importances()

		self.mean_tpr = mean_tpr
		self.mean_fpr = mean_fpr
		print('=====================\n\n')
	
		self.metrics_list = metrics_list
		#print("\nMetrics: ", metrics_list)
	
	
	
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
		elif self.model_type == 'RF':
			importances = self.model.feature_importances_
		
		
		importances_series = pd.Series(importances, index=self.feature_cols)
		importances_series.sort_values(inplace=True)
		print(importances_series)
		
		fig, ax = plt.subplots(figsize=(14, 10))
		
		importances_series.plot.bar()
		plt.show()

		pdf_filename = self.out_dir + '/RF_importance_scores.' + self.score_print_name + '.' + self.model_type + '.' + self.Y_label

		if self.include_vcf_extracted_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')
