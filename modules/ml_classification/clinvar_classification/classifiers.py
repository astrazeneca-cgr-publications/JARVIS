import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV, LinearRegression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_curve, auc, mean_squared_error, mean_absolute_error, explained_variance_score, r2_score
from sklearn.model_selection import train_test_split
import sys
import os

sys.path.insert(1, os.path.join(sys.path[0], '../..'))
from custom_utils import create_out_dir
	


class Classifier:
	
	def __init__(self, Y_label, out_dir, regression=False, model_type='RandomForest', 			 include_vcf_extracted_features=False, base_score='gwrvis', use_only_base_score=False):
		self.out_dir = out_dir
		self.Y_label = Y_label
		self.regression = regression
		self.model_type = model_type
		self.include_vcf_extracted_features = include_vcf_extracted_features
		self.base_score = base_score
		self.use_only_base_score = use_only_base_score

		

		print('\n-----\nRegression:', self.regression)
		print('Model:', model_type, '\n-----\n')
		
		
	
	def split_into_train_test_sets(self, df):

	
		if self.use_only_base_score:
		
			df.dropna(inplace=True)
			self.feature_cols = [self.base_score]
			
		else:
			cols_to_drop = ['chr', 'start', 'end', 'genomic_class', self.Y_label, 'clinvar_annot']
	
			vcf_dependent_cols = ['common_variants', 'common_vs_all_variants_ratio', 'all_variants', 'mean_ac', 'mean_af', 'bin_1', 'bin_2', 'bin_3', 'bin_4', 'bin_5', 'bin_6']
			
			if not self.include_vcf_extracted_features:
				cols_to_drop.extend(vcf_dependent_cols)
			
			self.feature_cols = [x for x in df.columns.values if x not in cols_to_drop ]
			print('Features:', self.feature_cols)
			# BETA
			# Retaining only most important features (improves AUC only by +0.001)
			#self.feature_cols = ['gwrvis', 'gc_content', 'cpg', 'mut_rate', 'cpg_islands', 'H3K4me2']
			
			
		fig, ax = plt.subplots(figsize=(10,10))	
		pathogenic = df.loc[ df.clinvar_annot == 1, self.base_score].tolist()
		benign = df.loc[ df.clinvar_annot == 0, self.base_score].tolist()
		
		sns.distplot(pathogenic, hist=False, kde=True, label='pathogenic (' + str(len(pathogenic)) + ')')
		sns.distplot(benign, hist=False, kde=True, label='benign (' + str(len(benign)) + ')')
		plt.title(self.base_score)
		plt.close()
		fig.savefig('cadd_intergenic.density_plots.pdf')
				
		print(self.feature_cols)
		print(self.Y_label)
		
		
		X = df[self.feature_cols]
		y = df[self.Y_label].astype(int).values

		print(X)
		print(y)
		
		#self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)

		# with CV
		self.X_train = X
		self.X_test = X 
		self.y_train = y 
		self.y_test = y
		
		print(self.X_train)
		print(self.y_train)
		
		
	def run_classifier_model(self):

		rf_params = dict(n_estimators=100, max_depth=2, random_state=0)
	
		if self.regression:
			# linear regression
			if self.model_type == 'RandomForest':
				model = RandomForestRegressor(**rf_params)
				
			elif self.model_type == 'Logistic':
				model = LinearRegression()
			
			model.fit(self.X_train, self.y_train)
			
			self.y_pred = model.predict(self.X_test)
			rmse = mean_squared_error(self.y_test, self.y_pred)
			print('RMSE:', rmse)
			mae = mean_absolute_error(self.y_test, self.y_pred)
			print('Mean Absolute Error:', mae)
			evs = explained_variance_score(self.y_test, self.y_pred)
			print('Explained variance score:', evs)
			
			# r2 represents the proportion of variance (of y) that has been explained by the independent variables in the model
			r2 = r2_score(self.y_test, self.y_pred)
			print('R^2:', r2)
			
		else:
			# classification
			if self.model_type == 'RandomForest':
				model = RandomForestClassifier(**rf_params)	
				#add cross-validation
				print(np.mean(cross_val_score(model, self.X_train, self.y_train, cv=5)))
				sys.exit()
				
			elif self.model_type == 'Logistic':
				print('> Running logistic regression...')
				#model = LogisticRegression(C=1e9, solver='lbfgs', max_iter=10000)
				model = LogisticRegressionCV(cv=5, solver='lbfgs', max_iter=10000)
			
			model.fit(self.X_train, self.y_train)
			
			self.pred_proba = model.predict_proba(self.X_test)[:, 1]
		
		return model
	
		

	def plot_roc_curve(self, model):

		fpr, tpr, thresholds = roc_curve(self.y_test, self.pred_proba)
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
		

		pdf_filename = self.out_dir + '/' + self.model_type + '_ROC_' + self.base_score + \
						'.AUC_' + str(roc_auc)
						
		if self.use_only_base_score:
			pdf_filename += '.use_only_base_score'
		if self.include_vcf_extracted_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')
		
		
	
	def get_feature_importances(self, model):
	
		print("\n> Feature importances:")
		if self.model_type == 'Logistic':
			importances = model.coef_.reshape(-1, 1)[:, 0]
		elif self.model_type == 'RandomForest':
			importances = model.feature_importances_
		
		
		importances_series = pd.Series(importances, index=self.feature_cols)
		importances_series.sort_values(inplace=True)
		print(importances_series)
		
		fig, ax = plt.subplots(figsize=(14, 10))
		
		importances_series.plot.bar()
		plt.show()

		pdf_filename = self.out_dir + '/RF_importance_scores.' + self.base_score + '.' + self.model_type

		if self.use_only_base_score:
			pdf_filename += '.use_only_base_score'
		if self.include_vcf_extracted_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')
		
		
	def train_and_predict(self):
	
		model = self.run_classifier_model()
		
		if not self.regression:
			self.plot_roc_curve(model)

		self.get_feature_importances(model)