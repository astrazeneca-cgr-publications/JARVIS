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
from prepare_data import JarvisDataPreprocessing

import sys, os 
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
sys.path.insert(1, os.path.join(sys.path[0], '..')) 
sys.path.insert(1, os.path.join(sys.path[0], '.')) 
import custom_utils



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
	
		print('indexes:', len(indexes))		
		for key, _ in data_dict.items():
			data_dict[key] = data_dict[key][indexes]
			print(key, ':', len(data_dict[key]))

		return data_dict


	def filter_data_by_genomic_class(self, data_dict, genomic_classes=['intergenic']):
		
		print('Retaining only variants for class:', ', '.join(genomic_classes)) 

		all_class_indexes = np.array([])
		for cur_class in genomic_classes:
			cur_class_indexes = np.argwhere(data_dict['genomic_classes'] == cur_class)
			print(cur_class, ':', len(cur_class_indexes))
			all_class_indexes = np.append(all_class_indexes, cur_class_indexes)

		all_class_indexes = all_class_indexes.astype(int).tolist()
		print(all_class_indexes[:10])
		print(len(all_class_indexes))
		
		filtered_data_dict = self.subset_data_dict_by_index(data_dict, all_class_indexes)
		print('filtered unique genomic_classes:', np.unique(filtered_data_dict['genomic_classes']))
		sys.exit()

		return filtered_data_dict



	def compile_and_train_model(self, model, data_dict):
		
		learning_rate = 0.001 # default: 0.001
		adam_optimizer = tf.keras.optimizers.Adam(learning_rate=learning_rate)

		model.compile(loss='binary_crossentropy', optimizer=adam_optimizer, metrics=[tf.keras.metrics.CosineSimilarity(axis=1)]) #try 'binary_accuracy' for metrics


		checkpoint_name = 'gwrvis_best_model.hdf5'
		checkpointer = ModelCheckpoint(checkpoint_name, monitor='val_loss', verbose=1, save_best_only=True, mode='auto')
		earlystopper = EarlyStopping(monitor='val_loss', patience=10, verbose=1)

		if self.include_vcf_features:
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
	include_vcf_features = False

	jarvis_trainer = JarvisTraining(config_file, include_vcf_features)

	


	# Read train, validation, test data
	data_dict = jarvis_trainer.read_data()

	# Print input data shapes for sanity check check for consistency
	jarvis_trainer.inspect_input_data(data_dict)

	jarvis_trainer.filter_data_by_genomic_class(data_dict, ['utr', 'intergenic'])
	sys.exit()


	## NEED TO ADD MODULE TO FILTER BY GENOMIC CLASS

	## NEED TO ADD module for cross validation
		## Function for selecting data based on indexes


	#model = nn_models.cnn_1_conv_2_fcc(win_len) 
	model = nn_models.cnn_2_conv_2_fcc(win_len) 
	#model = nn_models.cnn_rnn_1_conv_1_lstm(win_len)
	
	#model = functional_nn_models.funcapi_cnn_1_conv_2_fcc(win_len)
	print(model.summary())



	model, history = compile_and_train_model(model, train_dict, validation_dict, include_ext_feat=include_ext_feat)

	plot_history(history)

	test_and_evaluate_model(model, test_dict, include_ext_feat=include_ext_feat)
