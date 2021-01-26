from Model import Model
import tensorflow as tf
from utils import *
from tensorflow.keras.models import load_model
from keras.models import Sequential
from keras.layers import Dense, Activation, Embedding, Flatten
from keras.layers import LSTM, Input, GRU
from keras.layers.normalization import BatchNormalization
from keras import optimizers
from sklearn.decomposition import PCA
import numpy as np
import pickle

def get_shalla_model():
	return GRUModel('../lib/glove.6B.50d-char.txt', 50, pca_embedding_dim=16, maxlen=50, gru_size=16, batch_size=8192, lr=0.005, decay=0.001, hidden_size=8, epochs=3, model_dir="../lib/shalla/")

def get_ycsb_model():
	return GRUModel('../lib/glove.6B.50d-char.txt', 23, pca_embedding_dim=None, maxlen=23, lossfc="mse", batch_size=100000, dense_only=True, lr=0.01, hidden_size=24, epochs=10, model_dir="../lib/ycsb/")

class GRUModel(Model):
	def __init__(self, embeddings_path, embedding_dim, lr=0.001, maxlen=50, pca_embedding_dim=None, lossfc="binary_crossentropy", batch_size=1024, gru_size=16, 
	hidden_size=None, second_gru_size=None, decay=0.001, epochs=10, lstm=False, dense_only=False, model_dir="../lib/shalla/"):
		self.embeddings_path = embeddings_path
		self.embedding_dim = embedding_dim
		self.lr = lr
		self.maxlen = maxlen
		self.pca_embedding_dim = pca_embedding_dim
		self.model = None
		self.batch_size = batch_size
		self.gru_size = gru_size
		self.hidden_size = hidden_size
		self.second_gru_size = second_gru_size
		self.decay = decay
		self.epochs = epochs
		self.lstm = lstm
		self.lossfc = lossfc
		self.dense_only = dense_only
		self.model_dir = model_dir

	def initmodel(self):
		num_chars = len(self.char_indices)
		embedding_vectors = {}
		with open(self.embeddings_path, 'r') as f:
		    for line in f:
		        line_split = line.strip().split(" ")
		        vec = np.array(line_split[1:self.maxlen+1], dtype=float)
		        char = line_split[0]
		        embedding_vectors[char] = vec

		embedding_matrix = np.zeros((num_chars + 1, self.embedding_dim))
		for char, i in self.char_indices.items():
		    embedding_vector = embedding_vectors.get(char)
		    assert(embedding_vector is not None)
		    embedding_matrix[i] = embedding_vector

		print(embedding_matrix.shape)

		if self.pca_embedding_dim:
		    pca = PCA(n_components=self.pca_embedding_dim)
		    pca.fit(embedding_matrix[1:])
		    embedding_matrix_pca = np.array(pca.transform(embedding_matrix[1:]))
		    embedding_matrix_pca = np.insert(embedding_matrix_pca, 0, 0, axis=0)
		    print("PCA matrix created")
		
		if not self.lstm:
			rnn_layer = GRU(self.gru_size, return_sequences=False if not self.second_gru_size else True)
		else:
			rnn_layer = LSTM(self.gru_size, return_sequences=False if not self.second_gru_size else True)


		prelayers = [
			Embedding(num_chars + 1, self.embedding_dim if not self.pca_embedding_dim else self.pca_embedding_dim, input_length=self.maxlen,
    weights=[embedding_matrix] if not self.pca_embedding_dim else [embedding_matrix_pca]),
		    # Embedding(num_chars + 1, self.embedding_dim if not self.pca_embedding_dim else self.pca_embedding_dim, input_length=self.maxlen),
		    rnn_layer
		]

		if self.second_gru_size:
			prelayers.append(GRU(self.second_gru_size))

		if self.hidden_size:
			prelayers.append(Dense(self.hidden_size)),
			prelayers.append(Activation('relu'))

		postlayers = [
			Dense(1),
			Activation('sigmoid'),
		]

		if not self.dense_only:
			layers = prelayers + postlayers
		else:
			layers = [
		    Embedding(num_chars + 1, self.embedding_dim if not self.pca_embedding_dim else self.pca_embedding_dim, input_length=self.maxlen,
    weights=[embedding_matrix] if not self.pca_embedding_dim else [embedding_matrix_pca]),
			# Embedding(num_chars + 1, self.embedding_dim if not self.pca_embedding_dim else self.pca_embedding_dim, input_length=self.maxlen),
		    Flatten(),
			Dense(self.hidden_size*16),
			BatchNormalization(),
		    Activation('relu'),
			Dense(self.hidden_size*8),
			BatchNormalization(),
		    Activation('relu'),
		    Dense(self.hidden_size*4),
			BatchNormalization(),
		    Activation('relu'),
		    Dense(self.hidden_size*2),
			BatchNormalization(),
		    Activation('relu'),
		    Dense(self.hidden_size),
			BatchNormalization(),
		    Activation('sigmoid'),
		    Dense(1),
		    Activation('sigmoid'),
		]
		self.model = Sequential(layers)
		optimizer = optimizers.Adam(lr=self.lr, decay=self.decay)
		self.model.compile(optimizer=optimizer, loss=self.lossfc, metrics=['accuracy'])

	def fit(self, text_X, text_y):
		X, y, self.char_indices, self.indices_char = vectorize_dataset(text_X, text_y, self.maxlen)
		f = open(self.model_dir + "char_indices.pkl", 'wb')
		pickle.dump(self.char_indices, f)
		f.close()
		self.initmodel()
		self.model.fit(X, y, batch_size=self.batch_size, epochs=self.epochs, verbose=2)
		self.model.save(self.model_dir + "model.h5")

	def predict(self, text_x):
		x = np.zeros((1, self.maxlen), dtype=np.int)
		offset = max(self.maxlen - len(text_x), 0)
		for t, char in enumerate(text_x):
		    if t >= self.maxlen:
		        break
		    x[0, t + offset] = self.char_indices[char]
		pred = self.model.predict(x)
		return pred[0][0]

	# Like predict, but you pass in an array of URLs, and it is all
	# vectorized in one step, making it more efficient
	def predicts(self, text_X):
		X = np.zeros((len(text_X), self.maxlen), dtype=np.int)
		for i in range(len(text_X)):
			offset = max(self.maxlen - len(text_X[i]), 0)
			for t, char in enumerate(text_X[i]):
			    if t >= self.maxlen:
			        break
			    X[i, t + offset] = self.char_indices[char]
		preds = [pred[0] for pred in self.model.predict(X, batch_size=8192)]
		return preds
		
	def loadmodel(self):
		f = open(self.model_dir + "char_indices.pkl", 'rb')
		self.char_indices = pickle.load(f)
		f.close()
		self.initmodel()
		self.model = load_model(self.model_dir + "model.h5")

	def getFeatures(self, text_X):
		X = np.zeros((len(text_X), self.maxlen), dtype=np.int)
		for i, url in enumerate(text_X):
			offset = max(self.maxlen - len(url), 0)
			for t, char in enumerate(url):
				if t >= self.maxlen:
					break
				X[i, t + offset] = self.char_indices[char]
		return X
