'''A neural network used to predict future internal peptide coordinates.'''

import keras
import matplotlib.pyplot as plt
from keras.layers import Dense, Input
from keras.models import Sequential

class Predictor():
    '''A basic feed-forward neural net to predict future internal peptide coordinates'''

    def __init__(self, verbose=False):
        super().__init__()
        model = _construct_model(input_dim=8, output_dim=5, hidden_layers=[5, 5, 5])
        self.model = _compile_model(model)
        self.verbose = verbose

    def run(self, x_train, y_train):
        if self.verbose:
            print('Training Model...')
        history = self.model.fit(x_train, y_train, epochs=100, batch_size=64)
        # history = self.model.fit(x_train, y_train, validation_data=(x_test, y_test), 
        #         epochs=100, batch_size=64)

        # plt.plot(history.history['acc'])
        # plt.plot(history.history['val_acc'])
        # plt.title('Model accuracy')
        # plt.ylabel('Accuracy')
        # plt.xlabel('Epoch')
        # plt.legend(['Train', 'Test'], loc='upper left')
        # plt.show()

        # plt.plot(history.history['loss'])
        # plt.plot(history.history['val_loss']) 
        # plt.title('Model loss') 
        # plt.ylabel('Loss') 
        # plt.xlabel('Epoch') 
        # plt.legend(['Train', 'Test'], loc='upper left') 
        # plt.show()

    def _compile_model(self, model):
        '''Specify the loss function and optimizer for our model.
        
        ### Loss Functions:
        - `SGD` (stochastic gradient descent)
        - `RMSProp` (root mean squared prop, see Hinton 2012)
        - `Adam` (RMSProp with momentum)
        - `Adadelta`
        - `Adagrad`
        - `NAdam` (Adam with nesterov momentum)
        - `FTRL`

        ### Regression Losses:
        - `mean_squared_error` (computes mean squares of errors)
        - `mean_absolute_error` (computes mean absolute percentage error)

        ### Metrics:
        - `root_mean_squared_error`
        - `mean_squared_error`
        '''
        if self.verbose:
            print('Compiling Model...')
        model.compile(
            loss='mean_squared_error',
            optimizer='sgd',
            metrics=['root_mean_squared_error']
        )
        return model

    def _construct_model(self, input_dim, output_dim, hidden_layers, activation_function='relu'):
        '''Constructs and returns the model.
        
        Parameters
        ----------
        input_dim - the number of dimensions of the input layer
        output_dim - the number of dimensions of the output layer
        hidden_layers - a list of the dimension of each hidden layer where the 
            length of list represents the number of hidden layers.
        activation_function - the activation function to use between layers (default RELU)
        '''
        if self.verbose:
            print('Constructing Model...')

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
