import matplotlib
matplotlib.use('agg') 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle

from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_auc_score

import logging 
logging.getLogger('tensorflow').disabled = True
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Convolution1D, MaxPooling1D
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
import nn_models
import functional_nn_models

import sys, os 
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
sys.path.insert(1, os.path.join(sys.path[0], '..')) 
import custom_utils


def init_dirs(config_file):
	out_dir = custom_utils.create_out_dir(config_file)
	out_dir = '../' + out_dir	
	ml_data_dir = out_dir + '/ml_data'

	return out_dir, ml_data_dir


def read_data(random_seqs=False):
	print('> Loading data (top ' + str(100 * float(top_ratio)) + ' %) ...')
	top_ratio_str = '.top_' + str(top_ratio)

	if random_seqs:
		print('[!] Reading data from random sequences...')
		top_ratio_str += '.random'

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
	


def compile_and_train_model(model, train_dict, validation_dict, include_ext_feat=False):
	
	learning_rate = 0.001 # default: 0.001
	adam_optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

	if regression:
		model.compile(loss='mean_absolute_error', optimizer=adam_optimizer, metrics=['mean_absolute_error'])
	else:
		model.compile(loss='binary_crossentropy', optimizer=adam_optimizer, metrics=[tf.keras.metrics.CosineSimilarity(axis=1)]) #try 'binary_accuracy' for metrics


	checkpoint_name = 'gwrvis_best_model.hdf5'
	checkpointer = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
	earlystopper = EarlyStopping(monitor='val_loss', patience=10, verbose=1)

	if include_ext_feat:
		history = model.fit([train_dict['seqs'], train_dict['X_other']], [train_dict[y]], batch_size=nn_models.batch_size, epochs=1000, 
			  shuffle=True,
			  validation_data=([validation_dict['seqs'], validation_dict['X_other']], [validation_dict[y]]), 
			  callbacks=[checkpointer,earlystopper])
	else:
		history = model.fit(train_dict['seqs'], train_dict[y], batch_size=nn_models.batch_size, epochs=1000, 
			  shuffle=True,
			  validation_data=(validation_dict['seqs'], validation_dict[y]), 
			  callbacks=[checkpointer,earlystopper])

	return model, history
	

def plot_history(history):

	fig, ax = plt.subplots(figsize=(15, 15))

	plt.plot(history.history['loss'], label='train')
	plt.plot(history.history['val_loss'], label='validation')
	plt.title('Model loss')
	plt.ylabel('loss')
	plt.xlabel('epoch')
	plt.legend()
	
	fig.savefig(out_dir + '/model_loss_history.pdf', bbox_inches='tight')
	
	

def get_metrics(test_flat, preds_flat):

	roc_auc = roc_auc_score(test_flat, preds_flat)
	

	accuracy = accuracy_score(test_flat, preds_flat)
	confus_mat =  confusion_matrix(test_flat, preds_flat) 
	print('> Confusion matrix:\n', confusion_matrix(test_flat, preds_flat))
		
	TN, FP, FN, TP = confus_mat.ravel()	
	print('TN:', TN, '\nFP:', FP, '\nFN:', FN, '\nTP:', TP)


	sensitivity = TP / (TP + FN)
	precision = TP / (TP + FP)
	specificity = TN / (TN + FP)

	print('> ROC AUC:', roc_auc)
	print('> Accuracy:', accuracy_score(test_flat, preds_flat))

	print('\n> Sensitivity:', sensitivity)
	print('> Precision:', precision)
	print('> Specificity:', specificity)
	
	metrics = {'roc_auc': roc_auc, 'accuracy': accuracy, 'sensitivity': sensitivity, 'precision': precision, 'specificity': specificity, 'confusion_matrix': confus_mat}

	return metrics


		  
def test_and_evaluate_model(model, test_dict, include_ext_feat=False):

	if include_ext_feat:
		print(">> Test: Including external features...")
		test_results = model.evaluate([test_dict['seqs'], test_dict['X_other']], [test_dict[y]]) 
		preds = model.predict([test_dict['seqs'], test_dict['X_other']])
	else:
		test_results = model.evaluate(test_dict['seqs'], test_dict[y]) 
		preds = model.predict(test_dict['seqs'])
	print(test_results)
	
	decision_thres = 0.5 # for classification
	if regression:
		decision_thres = 0 # 0 is the natural border between tolerant and intolerant gwRVIS values

	preds[preds >= decision_thres] = 1
	preds[preds < decision_thres] = 0
	
	preds_flat = np.argmax(preds, axis=1)
	test_flat = np.argmax(test_dict[y], axis=1)
	
	metrics = get_metrics(test_flat, preds_flat)


def compute_salient_bases_in_seq(model, seq):

	input_tensors = [model.input]
	gradients = model.optimizer.get_gradients(model.output[0][1], model.input)
	compute_gradients = K.function(inputs = input_tensors, outputs = gradients)
	x_value = np.expand_dims(seq, axis=0)
	gradients = compute_gradients([x_value])[0][0]
	sal = np.clip(np.sum(np.multiply(gradients, seq), axis=1), a_min=0, a_max=None)   

	return sal
	

def get_avg_saliency_scores(model, sequences):
	
	all_sal = []

	# TODO: there's no point aggregating scores from all sequences at each base index.
	# Maybe, for each sequence extract a contiguous substring that has the highest saliency scores
	# (allowing 1 or 2 'gaps' with lower saliency scores) 
	# and then do multiple sequence alighment between the corresponding substrings to see
	# if there is any particular motif.

	for seq in sequences:
		sal = compute_salient_bases(model, input_features[sequence_index])


	fig, ax = plt.subplots(figsize=(15, 15))

	barlist = plt.bar(np.arange(len(sal)), sal)
	plt.xlabel('Bases')
	plt.ylabel('Magnitude of saliency values')
	plt.xticks(np.arange(len(sal)), list(sequences[sequence_index]));

	plt.title('Avg. Saliency map for bases across all sequences')

	fig.savefig(out_dir + '/avg_saliency_per.pdf', bbox_inches='tight')

	


	
if __name__ == '__main__':

	config_file = sys.argv[1]
	top_ratio = sys.argv[2] 	#default: 0.01 -- look at top 1% of intolerant/tolerant windows
	random_seqs = bool(int(sys.argv[3])) # 1 for True or 0 for False

	# Parameters
	regression=False
	include_ext_feat = True


	if regression:
		y = 'gwrvis'   # continuous value
	else:
		y = 'y'  # 1/0 annotation

	config_params = custom_utils.get_config_params(config_file)
	win_len = config_params['win_len']
	

	# Init out dir structure
	out_dir, ml_data_dir = init_dirs(config_file)

	# Read train, validation, test data
	train_dict, validation_dict, test_dict = read_data(random_seqs=random_seqs)
	print(validation_dict.keys())
	print('X:', validation_dict['X'].shape)
	print('y:', validation_dict['y'].shape)
	print('seqs:', validation_dict['seqs'].shape)
	print('gwrvis:', validation_dict['gwrvis'].shape)
	print('X_other:', validation_dict['X_other'].shape)

	# Print input data shapes to check for consistency
	inspect_input_data(train_dict, validation_dict, test_dict)

	#model = nn_models.cnn_1_conv_2_fcc(win_len, regression=regression) 
	model = nn_models.cnn_2_conv_2_fcc(win_len, regression=regression) 
	#model = nn_models.cnn_rnn_1_conv_1_lstm(win_len, regression=regression)
	
	#model = functional_nn_models.funcapi_cnn_1_conv_2_fcc(win_len, regression=regression)
	print(model.summary())
	#sys.exit()



	model, history = compile_and_train_model(model, train_dict, validation_dict, include_ext_feat=include_ext_feat)

	plot_history(history)

	test_and_evaluate_model(model, test_dict, include_ext_feat=include_ext_feat)
