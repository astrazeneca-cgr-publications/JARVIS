import matplotlib
matplotlib.use('agg') 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
import pickle
import uuid

from scipy import interp
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score
from sklearn.metrics import roc_curve, auc


import logging 
logging.getLogger('tensorflow').disabled = True
import sys, os 
#os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Convolution1D, MaxPooling1D
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam


from sklearn.model_selection import StratifiedKFold, cross_val_score
from imblearn.under_sampling import RandomUnderSampler
import nn_models
from func_api_nn_models import *
from prepare_data import JarvisDataPreprocessing



os.environ['KMP_WARNINGS'] = 'off'
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
#sys.path.insert(1, os.path.join(sys.path[0], '..')) 
sys.path.insert(1, os.path.join(sys.path[0], '.')) 
import custom_utils

	
	

def train_with_cv(X, input_features='structured', seqs=None, X_unscaled=True):
	""" 
		Arg 'input_features' may take 1 of 3 possible values:
		- stuctured: using only structured features as input
		- sequence: using only sequence features as input
		- both: using both structured and sequence features as inputs
	"""
	
	"""
	# Init output dir structure
	config_params = custom_utils.get_config_params(config_file)
	
	self.out_dir = custom_utils.create_out_dir(self.config_file)
	self.ml_data_dir = self.out_dir + '/ml_data'
	self.clinvar_feature_table_dir = self.ml_data_dir + '/clinvar_feature_tables'
	self.out_models_dir = self.ml_data_dir + '/models'
	"""
	
	if X_unscaled:
		pre_trained_model_path = '/projects/cgr/users/kclc950/JARVIS/out/topmed-NEW_X_unscaled_ClinVar_pathogenic-denovodb_benign-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/models/intergenic_utr_lincrna_ucne_vista'

		print('\nX: unscaled')
		
	else:
		pre_trained_model_path = '/projects/cgr/users/kclc950/JARVIS/out/topmed-NEW_ClinVar_pathogenic-denovodb_benign-winlen_3000.MAF_0.001.varType_snv.Pop_SNV_only-FILTERED/ml_data/models/intergenic_utr_lincrna_ucne_vista'

		print('\nX: normalised')

		
	pre_trained_model_file = pre_trained_model_path + '/JARVIS-' + input_features + '.model'
	
	pre_trained_model = load_model(pre_trained_model_file)
	print("Loaded pre-trained model:", pre_trained_model_file)

	
	
	# Test data
	print('- X features:', X.shape)
	print('- seqs:', seqs.shape)


	if input_features == 'structured':
		test_inputs = X

	elif input_features == 'sequences':
		test_inputs = seqs

	elif input_features == 'both':
		test_inputs = [X, seqs]
	



	
	# Get prediction probabilities per class 				
	probas_ = pre_trained_model.predict(test_inputs)  # for Keras functional API
	
	print("probas_:", probas_.shape)
	
	
	pred_scores = probas_[:, 1]
	print(pred_scores[:40])
	print(pred_scores[-40:])
	
	
	return pred_scores
	
	