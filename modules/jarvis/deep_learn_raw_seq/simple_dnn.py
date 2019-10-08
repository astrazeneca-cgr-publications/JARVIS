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





# ***************************************************************************
# ******			FEED FORWARD DNN 			*****
# ***************************************************************************
def create_feedf_dnn(input_dim, nn_arch=[32, 32]):

	print("\n\n> Creating feed-forward DNN model...")
	model = Sequential()

	layer_idx = 0
	for layer_size in nn_arch:
		if layer_idx == 0:
			model.add(Dense(nn_arch[layer_idx], input_dim=input_dim, activation='relu'))
		else:
			model.add(Dense(nn_arch[layer_idx], activation='relu'))
		layer_idx += 1

	model.add(Dense(2, activation='softmax'))	
	model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
	#model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[tf.keras.metrics.CosineSimilarity(axis=1)])
	print(model.summary())

	return model



def train_feedf_dnn(model, data_dict, include_vcf_features, verbose=1):

	checkpoint_name = 'jarvis_best_feedf_model.hdf5'
	checkpointer = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose=verbose, save_best_only=True, mode='auto')
	earlystopper = EarlyStopping(monitor='val_loss', patience=20, verbose=verbose)


	if include_vcf_features:
		X = np.concatenate(data_dict['X'], data_dict['vcf_features'], axis=1)
		print('X with vcf features:', X.shape)

	batch_size = 64 #default values: 64 (intergenic), 64 (utr)
	history = model.fit(data_dict['X'], data_dict['y'], batch_size=batch_size, epochs=100, 
		  shuffle=True,
		  validation_split=0.1, 
		  callbacks=[checkpointer,earlystopper],
		  verbose=verbose)

	return history



