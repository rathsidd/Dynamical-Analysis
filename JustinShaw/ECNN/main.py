"""Run this file to start the model"""

import pywt

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
        self.data = range(2000)
        batch_size = 4
        assert((batch_size & (batch_size - 1) == 0) and batch_size != 0)
        self.batch_size = batch_size
    
    
    def run(self):
        '''Run this method to start the program.'''
        # get the first batch_size items and perfrom dwt on fine data
        fine_data = self._get_next_batch()

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
        coarse_batch_size = self.batch_size + 1
        for batch in self.batch_size:
            foo
        print(fine_batch)

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
            return pywt.dwt(data, 'haar') # TODO Is this the right transform?
        else:
            raise ValueError('Data should contain a power of two number of elements')

model = Model()
model.run()