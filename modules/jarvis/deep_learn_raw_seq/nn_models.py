import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten
from tensorflow.keras.layers import Convolution1D, MaxPooling1D
from tensorflow.keras.layers import LSTM
from tensorflow.keras.layers import Input, Dense, Bidirectional

## Global parameters for training
batch_size=1024



def cnn_1_conv_2_fcc(win_len, regression=False):

	# =================  Heuristic-based efficient architectures (win_len=3000) =================
	# - conv: Conv1D
	# - mp: MaxPooling1D
	# - dp: Dropout

	# >>> conv[32@50]_mp[5@15]_fcc[128]


	# CNN
	model = Sequential()
	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				padding="valid", strides=1, 
				filters=32, kernel_size=50)) #starting default: 30, 60 (ok)
	model.add(MaxPooling1D(strides=5, pool_size=50))

	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				padding="valid", strides=1, 
				filters=64, kernel_size=50)) #starting default: 30, 60 (ok)
	model.add(MaxPooling1D(strides=5, pool_size=50))

	#model.add(Convolution1D(activation="relu", 
	#			input_shape=(win_len, 4), 
	#			padding="valid", strides=1, 
	#			filters=32, kernel_size=40)) #starting default: 30, 60 (ok)
	#model.add(MaxPooling1D(strides=5, pool_size=40))


	# FCC
	model.add(Flatten())
	model.add(Dense(256, activation='relu'))
	#model.add(Dense(128, activation='relu'))
	#model.add(Dropout(0.2))

	if regression:
		model.add(Dense(1, kernel_initializer='normal', activation='linear'))
	else:
		model.add(Dense(2, activation='softmax'))

	return model



def cnn_2_conv_2_fcc(win_len, regression=False):
	
	# =================  Heuristic-based efficient architectures (win_len=3000) =================
	# - conv: Conv1D
	# - mp: MaxPooling1D
	# - dp: Dropout
	
	# >>> $ conv[32@50]_mp[5@15]_conv[64@25]_mp[5@15]_fcc[128]
	#
	#     > ROC AUC: 0.734375 > Accuracy: 0.7419354838709677  > Sensitivity: 0.5 > Precision: 0.9375 > Specificity: 0.96875 
	# - Train -- loss: 0.5287 - cosine_similarity: 0.7967
	# - Test -- loss: 0.5283 - cosine_similarity: 0.7979
	#
	#  <<< Legacy results >>>:
	#     ROC AUC: 0.76625 > Accuracy: 0.7709677419354839  > Sensitivity: 0.62 > Precision: 0.8691588785046729 > Specificity: 0.9125
	# - Train -- loss: 0.4509 - cosine_similarity: 0.8338
	# - Test -- loss: 0.4926 - cosine_similarity: 0.8122
	# [Comments]: slightly over-fitting	

	# >>> $ conv[32@50]_mp[5@15]_conv[64@25]_mp[5@15]_fcc[128]_dp[0.2]
	#
	#     > ROC AUC: 0.7372916666666666 > Accuracy: 0.7419354838709677  > Sensitivity: 0.5933333333333334 > Precision: 0.8240740740740741 > Specificity: 0.88125
	# - Train -- loss: 0.5550 - cosine_similarity: 0.7839 
	# - Test -- loss: 0.5120 - cosine_similarity: 0.8033
	# [Comments]: slight under-fitting
	
	# >>> $ conv[32@50]_mp[5@15]_conv[64@25]_mp[5@15]_fcc[128]_dp[0.2]_fcc[64]_dp[0.2]
	#
	#     > ROC AUC: 0.7556250000000001 > Accuracy: 0.7612903225806451  > Sensitivity: 0.58 > Precision: 0.8877551020408163 > Specificity: 0.93125
	# - Train -- loss: 0.5380 - cosine_similarity: 0.7932 
	# - Test -- loss: 0.5012 - cosine_similarity: 0.8072
	# [Comments]: slight under-fitting

	# Error >>> conv[32@50]_mp[5@15]_dp[0.2]_conv[64@25]_mp[5@15]_dp[0.2]_fcc[128]
	#
	#     ROC AUC: 0.7522916666666667 > Accuracy: 0.7580645161290323  > Sensitivity: 0.5733333333333334 > Precision: 0.8865979381443299 > Specificity: 0.93125
	# - Train -- loss: 0.5281 - cosine_similarity: 0.7961
	# - Test -- loss: 0.5134 - cosine_similarity: 0.8030 
	# [Comments]: no over-fitting; better specificity and precision (which is good to have low FP rate)
	# [Error]: with the dropout layers the network can't learn (apart from a legacy run only?)	

	# CNN
	model = Sequential()
	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				#kernel_initializer=tf.keras.initializers.he_uniform(),  	#or tf.contrib.layers.xavier_initializer(uniform=False)
				padding="valid", strides=1, 
				filters=32, kernel_size=50))
	model.add(MaxPooling1D(strides=5, pool_size=15))
	#model.add(Dropout(0.2))

	model.add(Convolution1D(activation="relu", 
				#kernel_initializer=tf.keras.initializers.he_uniform(),
				padding="valid", strides=1, 
				filters=64, kernel_size=25))

	model.add(MaxPooling1D(strides=5, pool_size=15))
	#model.add(Dropout(0.2))

	
	# FCC
	model.add(Flatten())
	model.add(Dense(128, activation='relu')) #, kernel_initializer=tf.keras.initializers.he_uniform()))
	#model.add(Dropout(0.2))
	#model.add(Dense(64, activation='relu'))
	#model.add(Dropout(0.2))

	if regression:
		model.add(Dense(1, kernel_initializer='normal', activation='linear'))
	else:
		model.add(Dense(2, activation='softmax'))

	return model

	"""
	# CNN
	model = Sequential()
	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				padding="valid", strides=1, 
				filters=128, kernel_size=50))
	model.add(MaxPooling1D(strides=15, pool_size=15))
	model.add(Dropout(0.2))

	model.add(Convolution1D(activation="relu", 
				padding="valid", strides=1, 
				filters=64, kernel_size=25))

	model.add(MaxPooling1D(strides=15, pool_size=15))
	model.add(Dropout(0.2))

	
	# FCC
	model.add(Flatten())
	model.add(Dense(128, activation='relu'))

	if regression:
		model.add(Dense(1, kernel_initializer='normal', activation='linear'))
	else:
		model.add(Dense(2, activation='softmax'))

	return model
	"""
		  

def cnn_rnn_1_conv_1_lstm(win_len, regression=False):
	
	forward_lstm = LSTM(units=320, return_sequences=True)
	brnn = Bidirectional(forward_lstm)
	
	# CNN
	model = Sequential()
	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				padding="valid", strides=1, 
				filters=16, kernel_size=30))
	model.add(MaxPooling1D(strides=15, pool_size=15))
	model.add(Dropout(0.2))

	# RNN
	model.add(brnn)
	model.add(Dropout(0.5))

	model.add(Flatten())
	model.add(Dense(128, activation='relu'))

	if regression:
		model.add(Dense(1, kernel_initializer='normal', activation='linear'))
	else:
		model.add(Dense(2, activation='softmax'))

	return model
		  
	
if __name__ == '__main__':
	
	model = cnn_1_conv_2_fcc(win_len=3000)
	print(model.summary())

	model = cnn_2_conv_2_fcc(win_len=3000)
	print(model.summary())

