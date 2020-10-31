"""Run this file to start the model"""

import pywt
from package.parser import Parser
from package.batch_manager import BatchManager
from package.predictor import Predictor

class Model:
    '''
    Error Correcting Neural Network
    -------------------------------
    The Error-Correcting Neural Network class manages the network 
    architecture to accelerate molecular dynamic simulations.
    '''

    def __init__(self):
        '''Initializes the Model class.'''

        # Gather the timeseries data as a list of tuples of sin/cos phi/psy angles
        self.data = Parser.get_data_from_files(verbose=True)
        self.batcher = BatchManager(self.data, verbose=True)
        self.predictors = []

    def run(self, verbose=False):
        '''Run this method to start the program.'''
        # Make a list of lists, where each inner list corresponds to one channel of the NN
        x_trains = []
        y_trains = []

        # Grab the next batch of data
        if verbose:
            print('Performing DWT on data and preparing for model training...')
        for fine_data, coarse_data, output_data in self.batcher.next():
            fine_dwts = []
            coarse_dwts = []
            outputs = []
            # invert the tuple of lists to a list of lists then do dwt
            for grouped_fine_data in self._to_inverted_list(fine_data):
                fine_dwts.append(self.dwt(grouped_fine_data))
            for grouped_coarse_data in self._to_inverted_list(coarse_data):
                coarse_dwts.append(self.dwt(grouped_coarse_data))
            for grouped_output_data in self._to_inverted_list(output_data):
                outputs.append(grouped_output_data)
            for i in range(len(fine_dwts)):
                x_trains.append(list(fine_dwts[i] + coarse_dwts[i]))
                y_trains.append(list(outputs[i]))
        
        # Construct the right number of predictor objects
        if verbose:
            print('Constructing predictor models...')
        for i in range(len(x_trains)):
            self.predictors[i] = Predictor(verbose=True)

        # Run each predictor object with the associated data
        if verbose:
            print('Starting to train models....')
        for predictor in self.predictors:
            predictor.run(x_trains[i], y_trains[i])
        
        # Done with program!
        if verbose:
            print('Done.')

    def _get_next_batch(self, start):
        '''Pulls the next batch from data and returns 
        
        Parameters
        ----------
        `start` (int) the index from which to start learning.

        Exceptions
        ----------
        Raises an IndexError if there are not enough elements in data for next batch.
        '''
        # check to make sure there is enough data to get the next batch
        if len(self.data) < self.batch_size * (self.batch_size + 1):
            raise IndexError('Not enough elements in data to take next batch.')
        
        # Pull batch_size elements from the end of list for fine-scale dwt
        fine_batch = self.data[start : start + self.batch_size]
        fine_dwt = self.dwt(fine_batch)
        print(fine_batch)

    def _to_inverted_list(self, data):
        '''When given a tuple of lists, it returns an inverted list of lists.'''
        result = []
        for i in range(len(data)):
            print(len(data[i]))
            for j in range(len(data[i])):
                if len(result) == j:
                    result.append([])
                result[j].append(data[i][j])
        return result

    def dwt(self, data):
        """
        This function performs a discrete wavelet transform on an array of 2^n
        values, returning the result.

        Exceptions
        ----------
        Raises a ValueError if the length of values is not a power of two.

        Returns
        -------
        Returns the result of a harr wavelet transform on the input values.
        """
        # Use bit manipulations to check for power of two, for more information
        # see here: https://stackoverflow.com/a/57025941
        n = len(data)
        if (n & (n - 1) == 0) and n != 0:
            cA, cD = pywt.dwt(data, 'haar') # TODO Is this the right transform?
            return cA + cD
        else:
            print(data)
            raise ValueError('Data should contain a power of two number of elements')

model = Model()
model.run(verbose=True)