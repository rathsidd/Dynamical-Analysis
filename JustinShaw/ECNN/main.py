"""Run this file to start the model"""

from package.parser import Parser
from package.batch_manager import BatchManager
from package.predictors import Predictors

class Model:
    '''
    Error Correcting Neural Network
    -------------------------------
    The Error-Correcting Neural Network class manages the network 
    architecture to accelerate molecular dynamic simulations.
    '''

    def __init__(self):
        structure = Parser.get_structure_from_files(verbose=True, test_mode=True)
        batcher = BatchManager(structure, verbose=True)
        x, y = batcher.get_clean_dataset()
        self.predictors = Predictors(x, y,  verbose=True)

    def run(self):
        '''This method starts the program.'''
        self.predictors.predict(test_mode=True)

if __name__ == "__main__":
    model = Model()
    model.run()