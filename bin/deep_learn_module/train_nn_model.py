import numpy as np
import pandas as pd
import pickle

from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score
import keras
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution1D, MaxPooling1D
from keras.callbacks import ModelCheckpoint, EarlyStopping
import nn_models

import sys, os 
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
sys.path.insert(1, os.path.join(sys.path[0], '..')) 
import custom_utils


def init_dirs(config_file):
	out_dir = custom_utils.create_out_dir(config_file)
	out_dir = '../' + out_dir	
	ml_data_dir = out_dir + '/ml_data'

	return out_dir, ml_data_dir


def read_data():
	print('> Loading data (top ' + str(100 * float(top_ratio)) + ' %) ...')
	top_ratio_str = '.top_' + str(top_ratio)

	print('Reading training data...')
	pkl_in = open(ml_data_dir + '/train' + top_ratio_str + '.pkl', 'rb')
	train_dict = pickle.load(pkl_in)
	pkl_in.close()

	print('Reading validation data...')
	pkl_in = open(ml_data_dir + '/validation' + top_ratio_str + '.pkl', 'rb')
	validation_dict = pickle.load(pkl_in)
	pkl_in.close()

	print('Reading test data...')
	pkl_in = open(ml_data_dir + '/test' + top_ratio_str + '.pkl', 'rb')
	test_dict = pickle.load(pkl_in)
	pkl_in.close()
	
	return train_dict, validation_dict, test_dict


def inspect_input_data(train_dict, validation_dict, test_dict):

	print('\n> Train:')
	print(train_dict['X'].shape)
	print(train_dict[y].shape)
	print(train_dict['seqs'].shape)

	print('\n> Validation:')
	print(validation_dict['X'].shape)
	print(validation_dict[y].shape)
	print(validation_dict['seqs'].shape)

	print('\n> Test:')
	print(test_dict['X'].shape)
	print(test_dict[y].shape)
	print(test_dict['seqs'].shape)
	


def compile_and_train_model(model):
	
	if regression:
		model.compile(loss='mean_absolute_error', optimizer='adam', metrics=['mean_absolute_error'])
	else:
		model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['cosine'])

	checkpoint_name = 'gwrvis_best_model.hdf5'
	checkpointer = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
	earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)

	model.fit(train_dict['seqs'], train_dict[y], batch_size=2048, epochs=100, 
		  shuffle=True,
		  validation_data=(validation_dict['seqs'], validation_dict[y]), 
		  callbacks=[checkpointer,earlystopper])

	return model
		  
		  
def test_and_evaluate_model(model):
	test_results = model.evaluate(test_dict['seqs'], test_dict[y]) 
	print(test_results)
	
	preds = model.predict(test_dict['seqs'])
	decision_thres = 0.5 # for classification
	if regression:
		decision_thres = 0 # 0 is the natural border between tolerant and intolerant gwRVIS values

	preds[preds >= decision_thres] = 1
	preds[preds < decision_thres] = 0
	
	preds_flat = preds.flatten()
	test_flat = test_dict[y].flatten()
	
	print(accuracy_score(test_flat, preds_flat))
	print(confusion_matrix(test_flat, preds_flat))
	
	roc_auc = roc_auc_score(test_flat, preds_flat)
	print('ROC AUC:', roc_auc)


	
if __name__ == '__main__':

	config_file = sys.argv[1]
	top_ratio = sys.argv[2] 	#default: 0.01 -- look at top 1% of intolerant/tolerant windows

	regression=True
	if regression:
		y = 'gwrvis'   # continuous value
	else:
		y = 'y'  # 1/0 annotation

	config_params = custom_utils.get_config_params(config_file)
	win_len = config_params['win_len']

	# init out dir structure
	out_dir, ml_data_dir = init_dirs(config_file)

	# read train, validation, test data
	train_dict, validation_dict, test_dict = read_data()

	# print input data shapes to check for consistency
	inspect_input_data(train_dict, validation_dict, test_dict)

	#model = nn_models.cnn_1_conv_2_fcc(regression=regression) 
	model = nn_models.cnn_rnn_1_conv_1_lstm(regression=regression)
	print(model.summary())

	model = compile_and_train_model(model)

	test_and_evaluate_model(model)
