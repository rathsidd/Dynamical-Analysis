'''A neural network used to predict future internal peptide coordinates.'''

import keras
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
from keras.layers import Dense, Input, LSTM
from keras.models import Sequential

class Predictors():
    '''A basic feed-forward neural net to predict future internal peptide coordinates'''

    def __init__(self, X, Y, verbose=False, test_mode=False):
        super().__init__()
        # Metaparameters
        self.learning_rate = 0.00001
        self.num_epochs = 1000
        self.size_batch = 1
        self.val_split = 0.2
        self.hidden_layer_list = [[2],[2]] # [model, err_model]
        # Class data
        self.show_dialog = verbose
        self.test_mode = test_mode
        self.X = tf.convert_to_tensor(X)
        self.Y = tf.convert_to_tensor(Y)
        self.models = []
        self.err_models = []
        self.histories = []
        if self.test_mode:
            model, err_model = self._construct_models(
                input_data=X[0][0],
                output_data=Y[0][0],
                hidden_layers=self.hidden_layer_list)
            self.models.append(model)
            self.err_models.append(err_model)
        else:
            for i, _ in enumerate(X):
                model, err_model = self._construct_models(
                    input_data=X[i][0],
                    output_data=Y[i][0],
                    hidden_layers=self.hidden_layer_list)
                self.models.append(model)
                self.err_models.append(err_model)
        self._compile_models()
    
    def predict(self):
        '''Fits the model with the given data.'''
        if self.test_mode:
            history = self.models[0].fit(
                self.X[0],
                self.Y[0],
                epochs=self.num_epochs,
                batch_size=self.size_batch,
                validation_split=self.val_split)
            err_history = self.err_models[0].fit(
                self.X[0],
                self.Y[0],
                epochs=self.num_epochs,
                batch_size=self.size_batch,
                validation_split=self.val_split)
            self.histories.append(err_history)
        else:
            for i, model in enumerate(self.models):
                history = model.fit(
                    self.X[i],
                    self.Y[i],
                    epochs=self.num_epochs,
                    batch_size=self.size_batch,
                    validation_split=self.val_split)
                self.histories.append(history)
        self._show_history()

    def _compile_model(self):
        '''Specify the loss function and optimizer for our model.'''
        opt = keras.optimizers.Adamax(learning_rate=self.learning_rate)
        err_opt = keras.optimizers.Adamax(learning_rate=self.learning_rate)
        for model in self.models:
            model.compile(loss='mean_absolute_error', optimizer=opt)
        for model in self.err_models:
            model.compile(loss='mean_absolute_error', optimizer=err_opt)

    def _construct_models(self, input_data, output_data, hidden_layers):
        '''Constructs a new model with the functional API.'''
        # Normalize input layer
        normalizer = keras.layers.experimental.preprocessing.Normalization()
        normalizer.adapt(input_data)
        inputs = Input(shape=input_data.shape)
        normalized_inputs = normalizer(inputs)
        model = self._construct_singular_model(normalized_inputs, hidden_layers[0], output_data.shape[0])
        err_model = self._construct_singular_model(model.outputs, hidden_layers[1], output_data.shape[0])
        return model, err_model
    
    def _construct_singular_model(self, inputs, hidden_layers, output_shape, act_func='relu'):
        '''Constructs a singular model and returns it'''
        # Add hidden layers to the main model
        prev_layer = inputs
        for layer in hidden_layers:
            hidden_layer = Dense(layer)(prev_layer)
            prev_layer = hidden_layer
        outputs = Dense(output_shape, activation=act_func)(prev_layer)
        return keras.Model(inputs=inputs, outputs=outputs)

    
    def _show_history(self):
        '''Takes an array of keras history objects and averages out the data in each, displaying the results.'''
        list_of_losses = []
        list_of_val_losses = []
        for h in self.histories:
            list_of_losses.append(h.history['loss'])
            list_of_val_losses.append(h.history['val_loss'])
        loss = np.average(list_of_losses, axis=0)
        val_loss = np.average(list_of_val_losses, axis=0)
        # Summarize history for loss
        plt.plot(loss)
        plt.plot(val_loss)
        plt.title('Residue-Averaged Model Loss')
        plt.ylabel('Average Loss')
        plt.xlabel('Epoch')
        plt.legend(['train', 'test'], loc='upper left')
        # plt.show()
