'''A neural network used to predict future internal peptide coordinates.'''

import keras
import numpy as np
import matplotlib.pyplot as plt
from keras.layers import Dense, Input
from keras.models import Sequential

class Predictors():
    '''A basic feed-forward neural net to predict future internal peptide coordinates'''

    def __init__(self, X, Y, verbose=False):
        super().__init__()
        self.show_dialog = verbose
        self.X = X
        self.Y = Y
        self.models = []
        self.histories = []
        for _ in range(len(X)):
            model = self._construct_model(input_dim=8, output_dim=5, hidden_layers=[7, 6])
            self._compile_model(model)
            self.models.append(model)
    
    def predict(self, test_mode=False):
        '''Fits the model with the given data.'''
        if test_mode:
            history = self.models[0].fit(self.X[0], self.Y[0], epochs=100, batch_size=10, validation_split=0.2)
            self.histories.append(history)
        else:
            for i, model in enumerate(self.models):
                history = model.fit(self.X[i], self.Y[i], epochs=100, batch_size=1, validation_split=0.2)
                self.histories.append(history)
        self._show_history(self.histories)

    def _compile_model(self, model):
        '''Specify the loss function and optimizer for our model.'''
        model.compile(
            loss='mean_absolute_error',
            optimizer=keras.optimizers.Adam(learning_rate=0.001),
        )
        return model

    def _construct_model(self, input_dim, output_dim, hidden_layers, activation_function='relu'):
        '''Returns a new model with the given attributes.
        
        Parameters
        ----------
        input_dim - the number of dimensions of the input layer
        output_dim - the number of dimensions of the output layer
        hidden_layers - a list of the dimension of each hidden layer where the 
            length of list represents the number of hidden layers.
        activation_function - the activation function to use between layers (default RELU)
        '''

        # Build a basic, sequential feed-forward NN model
        model = Sequential()        
        # Calculate the input layer size
        input_layer = Input(shape=(input_dim,))
        model.add(input_layer)

        # Add the hidden layers
        for layer in hidden_layers:
            hidden_layer = Dense(layer, activation=activation_function)
            model.add(hidden_layer)
        
        # Add the output layer to the model and return the model
        output_layer = Dense(output_dim)
        model.add(output_layer)
        return model
    
    def _show_history(self, histories):
        '''Takes an array of keras history objects and averages out the data in each, displaying the results.'''
        list_of_losses = []
        list_of_val_losses = []
        for h in histories:
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
        plt.show()
