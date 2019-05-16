import keras
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.convolutional import Convolution1D, MaxPooling1D
from keras.layers.recurrent import LSTM
from keras.layers import Bidirectional


def cnn_1_conv_2_fcc(win_len=3000, regression=False):
	
	model = Sequential()
	model.add(Convolution1D(activation="relu", 
				input_shape=(win_len, 4), 
				padding="valid", strides=1, 
				filters=16, kernel_size=30))
	model.add(MaxPooling1D(strides=15, pool_size=15))
	model.add(Dropout(0.2))

	model.add(Flatten())
	model.add(Dense(128, activation='relu'))

	if regression:
		model.add(Dense(1, kernel_initializer='normal', activation='linear'))
	else:
		model.add(Dense(2, activation='softmax'))

	return model



def cnn_2_conv_2_fcc(win_len=3000, regression=False):
	
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

	model.add(Flatten())
	model.add(Dense(128, activation='relu'))

	if regression:
		model.add(Dense(1, kernel_initializer='normal', activation='linear'))
	else:
		model.add(Dense(2, activation='softmax'))

	return model
		  

def cnn_rnn_1_conv_1_lstm(win_len=3000, regression=False):
	
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
	
	model = cnn_1_conv_2_fcc()
	print(model.summary())

	model = cnn_2_conv_2_fcc()
	print(model.summary())

