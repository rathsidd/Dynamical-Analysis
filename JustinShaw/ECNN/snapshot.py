'''This class defines the data structure that represents one frame in the PDB trajectory.'''

class Snapshot():
    '''
    ## Snapshot
    A data structure that represents a single frame of the MD trajectory, composed of
    the sines and cosines of the torsion angles of the original molecular structure.
    '''

    def __init__(self, model):
        '''Serializes a BioPython model into a Snapshot object'''
        self.