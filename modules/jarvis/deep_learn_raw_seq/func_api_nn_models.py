import logging 
logging.getLogger('tensorflow').disabled = True
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Dropout, Flatten
from tensorflow.keras.layers import LSTM, Bidirectional
from tensorflow.keras.layers import concatenate
import sys


def feedf_dnn(input_dim, nn_arch=[32,32]):

	features_input = Input(shape=(input_dim, ))

	layer_idx = 0
	for layer_size in nn_arch:
		if layer_idx == 0:
			x = Dense(nn_arch[layer_idx], activation='relu')(features_input)
			layer_idx += 1
		else:
			x = Dense(nn_arch[layer_idx], activation='relu')(x)
	
	output = Dense(2, activation='softmax')(x)
	
	model = Model(inputs=features_input, outputs=output)

	return model



def cnn2_fc2(win_len, num_features=4):

	seq_input = Input(shape=(win_len, 4), name='seq_input')

	# ---- seqs
	x = Conv1D(activation="relu", input_shape=(win_len, 4), padding="valid", strides=1, filters=128, kernel_size=11)(seq_input)
	x = MaxPooling1D(strides=4, pool_size=4)(x)
	x = Dropout(0.2)(x)

	x = Conv1D(activation="relu", padding="valid", strides=1, filters=64, kernel_size=11)(x)
	x = MaxPooling1D(strides=4, pool_size=4)(x)
	x = Dropout(0.2)(x)
	
	x = Flatten()(x)
	x = Dense(128, activation='relu')(x)
	x = Dense(32, activation='relu')(x)
	#x = Dense(16, activation='relu')(x)
	
	output = Dense(2, activation='softmax')(x)


	model = Model(inputs=seq_input, outputs=output)

	return model



def cnn2_concat_dnn_fc2(feat_input_dim, nn_arch=[32,32], win_len=3000, num_features=4):

	seq_input = Input(shape=(win_len, 4), name='seq_input')

	# ---- sequences as features
	x1 = Conv1D(activation="relu", input_shape=(win_len, 4), padding="valid", strides=1, filters=128, kernel_size=11)(seq_input)
	x1 = MaxPooling1D(strides=4, pool_size=4)(x1)
	x1 = Dropout(0.2)(x1)

	x1 = Conv1D(activation="relu", padding="valid", strides=1, filters=64, kernel_size=11)(x1)
	x1 = MaxPooling1D(strides=4, pool_size=4)(x1)
	x1 = Dropout(0.2)(x1)
	
	x1 = Flatten()(x1)
	x1 = Dense(128, activation='relu')(x1)
	x1 = Dense(32, activation='relu')(x1)
	seq_output = x1


	# ---- structured features
	feat_input = Input(shape=(feat_input_dim, ), name='feat_input')

	layer_idx = 0
	for layer_size in nn_arch:
		if layer_idx == 0:
			x2 = Dense(nn_arch[layer_idx], activation='relu')(feat_input)
			layer_idx += 1
		else:
			x2 = Dense(nn_arch[layer_idx], activation='relu')(x2)
	feat_output = x2



	# ---- concatenate 
	x = concatenate([seq_output, feat_output])
	x = Dense(64, activation='relu')(x)
	x = Dense(128, activation='relu')(x)

	
	output = Dense(2, activation='softmax')(x)


	model = Model(inputs=[feat_input, seq_input], outputs=output)


	return model



def cnn_2_conv_2_fcc(win_len, regression=False):
	
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
		  

def cnn_rnn_1_conv_1_lstm(win_len, regression=False):
	
	
	# CNN
	model = Sequential()
	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				padding="valid", strides=1, 
				filters=16, kernel_size=30))
	model.add(MaxPooling1D(strides=15, pool_size=15))
	model.add(Dropout(0.2))

	# RNN
	forward_lstm = LSTM(units=320, return_sequences=True)
	brnn = Bidirectional(forward_lstm)

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
	
	model = funcapi_cnn_1_conv_2_fcc(win_len=3000)
	print(model.summary())
	sys.exit()

	#model = cnn_2_conv_2_fcc()
	#print(model.summary())

