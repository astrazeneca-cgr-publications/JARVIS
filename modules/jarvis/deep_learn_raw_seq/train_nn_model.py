import matplotlib
matplotlib.use('agg') 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from scipy import interp
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score
from sklearn.metrics import roc_curve, auc


import logging 
logging.getLogger('tensorflow').disabled = True
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Convolution1D, MaxPooling1D
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import StratifiedKFold, cross_val_score
import nn_models
import func_api_nn_models
from prepare_data import JarvisDataPreprocessing

import sys, os 
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
#sys.path.insert(1, os.path.join(sys.path[0], '..')) 
sys.path.insert(1, os.path.join(sys.path[0], '.')) 
import custom_utils
from simple_dnn import create_feedf_dnn, train_feedf_dnn



class JarvisTraining:

	def __init__(self, config_file, include_vcf_features=False):
		
		self.config_file = config_file
		self.include_vcf_features = include_vcf_features		

		config_params = custom_utils.get_config_params(config_file)
		self.win_len = config_params['win_len']

		self.init_ouput_dirs()
		print(self.out_dir)


		self.data_dict_file = self.ml_data_dir + '/jarvis_data.pkl'
		if not os.path.exists(self.data_dict_file):
			print("Preparing training data - calling JarvisDataPreprocessing object...")
			data_preprocessor = JarvisDataPreprocessing(config_file)
			
			# Extract raw sequences from input variant windows and combine with original feature set         
			additional_features_df, filtered_onehot_seqs = data_preprocessor.compile_feature_table_incl_raw_seqs()        
			# Merge data, transform into form appropriate for DNN training and save into file
			self.data_dict_file = data_preprocessor.transform_and_save_data(additional_features_df, filtered_onehot_seqs)         
		print(self.data_dict_file)
		


	def init_ouput_dirs(self):

		# Init output dir structure
		self.out_dir = custom_utils.create_out_dir(self.config_file)
		self.ml_data_dir = self.out_dir + '/ml_data'


	def read_data(self):

		print('Reading training data...')
		pkl_in = open(self.data_dict_file, 'rb')
		data_dict = pickle.load(pkl_in)
		pkl_in.close()

		return data_dict


	def inspect_input_data(self, data_dict):

		print('\n> All data:')
		print('X:', data_dict['X'].shape)
		print('vcf_features:', data_dict['vcf_features'].shape)
		print('y:', data_dict['y'].shape)
		print('seqs:', data_dict['seqs'].shape)
		print('genomic_classes:', data_dict['genomic_classes'])
		print('genomic_classes:', type(data_dict['genomic_classes']))
		print('genomic_classes:', data_dict['genomic_classes'].shape)



	def subset_data_dict_by_index(self, data_dict, indexes):
	
		print('All indexes:', len(indexes))		
		for key, _ in data_dict.items():

			if key in ['X_cols', 'vcf_features_cols']:
				continue

			print('\n> Subsetting', key, ':', len(data_dict[key]))
			print(indexes[:10])
			print(data_dict[key][:10])

			data_dict[key] = data_dict[key][indexes]
			print('Filtered', key, ':', len(data_dict[key]))

		return data_dict


	def filter_data_by_genomic_class(self, data_dict, genomic_classes=None):
		
		print('Retaining only variants for class:', ', '.join(genomic_classes)) 

		all_class_indexes = np.array([])
		for cur_class in genomic_classes:
			cur_class_indexes = np.argwhere(data_dict['genomic_classes'] == cur_class)
			print(cur_class, ':', len(cur_class_indexes))
			all_class_indexes = np.append(all_class_indexes, cur_class_indexes)

		all_class_indexes = all_class_indexes.astype(int).tolist()
		
		filtered_data_dict = self.subset_data_dict_by_index(data_dict, all_class_indexes)
		print('filtered unique genomic_classes:', np.unique(filtered_data_dict['genomic_classes']))

		return filtered_data_dict	



	def include_vcf_features_for_prediction(self, data_dict):	
		X = np.concatenate(data_dict['X'], data_dict['vcf_features'], axis=1)
		
		return X



	def train_with_cv(self, data_dict, include_vcf_features, genomic_classes,
			  n_splits=5, batch_size=64, epochs=40, validation_split=0.1,
			  use_multiprocessing=True, verbose=0, input_features='structured'):
		""" 
			Arg 'input_features' may take 1 of 3 possible values:
			- stuctured: using only structured features as input
			- sequence: using only sequence features as input
			- both: using both structured and sequence features as inputs
		"""
		
		#if input_features != 'structured':
		#	verbose = 1


		if self.include_vcf_features:
			X = self.include_vcf_features_for_prediction(data_dict)
		else:
			X = data_dict['X'].copy()
		print('- X features:', X.shape)
		
		y = data_dict['y'].copy()
		
		if input_features != 'structured':
			seqs = data_dict['seqs'].copy()
		
	
		tprs = [] 
		aucs = [] 
		mean_fpr = np.linspace(0, 1, 100)  		
		fig, ax = plt.subplots(figsize=(10, 10))
	
		skf = StratifiedKFold(n_splits=n_splits, shuffle=True)

		i = 0
		for train_index, test_index in skf.split(X, np.argmax(y, axis=1)): #SK-fold requires non one-hot encoded y-data
			y_train, y_test = y[train_index], y[test_index]
			X_train, X_test = X[train_index], X[test_index]
			print('X_train:', X_train.shape)
			print('y_train:', y_train.shape)
			print('X_test:', X_test.shape)
			print('y_test:', y_test.shape)

			if input_features != 'structured':
				seqs_train, seqs_test = seqs[train_index], seqs[test_index]
				print('seqs_train:', seqs_train.shape)
				print('seqs_test:', seqs_test.shape)


			if input_features == 'structured':
				train_inputs = X_train
				test_inputs = X_test
			elif input_features == 'sequences':
				train_inputs = seqs_train
				test_inputs = seqs_test
			elif input_features == 'both':
				train_inputs = [X_train, seqs_train]
				test_inputs = [X_test, seqs_test]

			# -- Create new/clean model instance for each fold
			# > Keras functional API
			if input_features == 'structured':
				# @ Feed-forward DNN (for structured features input only)
				model = func_api_nn_models.feedf_dnn(filtered_data_dict['X'].shape[1], nn_arch=[64, 128, 256])	
			elif input_features == 'sequences':
				# @ CNN-CNN-FC-FC (for sequence features only)
				model = func_api_nn_models.cnn2_fc2(jarvis_trainer.win_len)
			model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])



			# ---- Callbacks ----
			checkpoint_name = 'jarvis_best_model_with_cv.hdf5'
			checkpointer = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose=verbose, save_best_only=True, mode='auto')
			patience = 10
			if input_features == 'structured':
				patience = 20
			earlystopper = EarlyStopping(monitor='val_loss', patience=patience, verbose=verbose)
			# -------------------


			history = model.fit(train_inputs, y_train, batch_size=batch_size, epochs=epochs, 
				  shuffle=True,
				  validation_split=validation_split, 
				  callbacks=[checkpointer, earlystopper], 
				  verbose=verbose)
			#self.plot_history(history, fold_id=(i+1))
			
			# Get prediction probabilities per class 				
			#probas_ = model.predict_proba(X_test)
			probas_ = model.predict(test_inputs)  # for Keras functional API

			# Compute ROC curve and area the curve 				
			fpr, tpr, thresholds = roc_curve(np.argmax(y_test, axis=1), probas_[:, 1])


			tprs.append(interp(mean_fpr, fpr, tpr)) 
			tprs[-1][0] = 0.0 
			roc_auc = round(auc(fpr, tpr), 3) 	
			print('Fold-', str(i+1), ' - AUC: ', roc_auc, '\n')
			aucs.append(roc_auc)
			plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
			
			i += 1

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
		plt.title('[' + ','.join(genomic_classes) + ']: ' + str(n_splits) + '-fold Cross-Validation ROC Curve')
		plt.legend(loc="lower right")
		#plt.show()
		#plt.close()
		
		
		pdf_filename = self.out_dir + '/Jarvis.' + input_features + '_features.' + '-'.join(genomic_classes) + '_ROC' + '.AUC_' + str(self.mean_auc)
		if include_vcf_features:
			pdf_filename += '.incl_vcf_features'
		pdf_filename += '.pdf'
		
		fig.savefig(pdf_filename, bbox_inches='tight')
		
		print('Mean AUC:', self.mean_auc)





	def plot_history(self, history, fold_id=0):

		fig, ax = plt.subplots(figsize=(15, 15))

		plt.plot(history.history['loss'], label='train')
		plt.plot(history.history['val_loss'], label='validation')
		plt.title('Model loss')
		plt.ylabel('loss')
		plt.xlabel('epoch')
		plt.legend()
		
		history_out_file = self.ml_data_dir + '/jarvis.model_loss_history.Fold-' + str(fold_id) + '.pdf'
		fig.savefig(history_out_file, bbox_inches='tight')
		
		print('Model training history saved into:', history_out_file)
		







if __name__ == '__main__':

	config_file = sys.argv[1]
	input_features = sys.argv[2]
	genomic_classes = sys.argv[3] # comma-separated
	genomic_classes = genomic_classes.split(',')	


	include_vcf_features = False




	jarvis_trainer = JarvisTraining(config_file, include_vcf_features)


	# Read train, validation, test data
	data_dict = jarvis_trainer.read_data()

	# Print input data shapes for sanity check check for consistency
	jarvis_trainer.inspect_input_data(data_dict)

	# Filter data by genomic class
	filtered_data_dict = jarvis_trainer.filter_data_by_genomic_class(data_dict, genomic_classes)
	print(filtered_data_dict)



	# Train using only structured features (i.e. without raw sequence data)	
	#model = create_feedf_dnn(filtered_data_dict['X'].shape[1], nn_arch=[32, 32]) 
	# ===== First test run - (only train set, without test set or CV) ==========
	#if False:
	#	history = train_feedf_dnn(model, filtered_data_dict, include_vcf_features)
	#	jarvis_trainer.plot_history(history)
	# ==========================================================================
	

	jarvis_trainer.train_with_cv(filtered_data_dict, include_vcf_features, genomic_classes, 
				     input_features=input_features, verbose=0)

