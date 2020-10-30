'''Generator that iterates and returns the next data for each batch.'''

import numpy as np

class BatchManager():

    def __init__(self, data, wavelet_size=4):
        '''Initialize the BatchManager class.
        
        This class manages a moving window over the given data based on the given wavelet size, and 
        '''
        super().__init__()
        self.data = data
        self.index = 0
        self.size = wavelet_size # this should be a power of two
        self.window_size = size + size * (size + 1)
    
    def next_batch(self):
        '''Returns a tuple of fine and coarse data.

        Returns
        -------
        The fine data is returned first as a list of (default 4) raw elements.
        The coarse data is returned second as a list of (default 4) averages of (default 5) elements. 
        '''
        try:
            if len(self.data) - self.index >= self.window_size:
                fine_data = self._get_fine_data()
                coarse_data = self._get_coarse_data()
                self.index += 1
                yield fine_data, coarse_data
        except:
            raise StopIteration

    def _get_fine_data(self):
        '''Returns a tuple of length `size` with elements that are themselves lists of sin/cos phi/psi angles.
        
        For example, this method might return:
        >>> ([model_1], [model_2], [model_3], [model_4])
        where each model is a list that contains data, for example `model_4` might be:
        >>> [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), sin(phi_2), cos(phi_2), ... , cos(psi_n)]
        '''
        result = ()
        models = self.data[self.index + self.window_size - self.size : self.index + self.window_size]
        for model in models:
            angles = self._get_angles_from_model()
            result += (angles)
        return result

    def _get_coarse_data(self):
        '''Returns a tuple of length `size` with elements that are themselves lists of the average of sin/cos phi/psi angles.'''
        result = ()
        models = self.data[self.index : self.index + self.window_size - self.size]
        for model in models:
            angles = self._get_angles_from_model()
            result += (angles)
        return self._take_average_of_tuple(result)

    def _get_angles_from_model(self, model):
        '''Returns a list of sin/cos phi/psi angles extracted from each model.'''
        angles = []
        for residue in model:
            sin_phi = residue.xtra['sin(phi)']
            cos_phi = residue.xtra['cos(phi)']
            sin_psi = residue.xtra['sin(psi)']
            cos_psi = residue.xtra['cos(psi)']
            # sin_omega = residue.xtra['sin(omega)']
            # cos_omega = residue.xtra['cos(omega)']
            # sin_chi2 = residue.xtra['sin(chi2)']
            # cos_chi2 = residue.xtra['cos(chi2)']
            angles.append(sin_phi) if sin_phi is not None else None
            angles.append(cos_phi) if cos_phi is not None else None
            angles.append(sin_psi) if sin_psi is not None else None
            angles.append(cos_psi) if cos_psi is not None else None
            # angles.append(sin_omega) if sin_omega is not None else None
            # angles.append(cos_omega) if cos_omega is not None else None
            # angles.append(sin_chi2) if sin_chi2 is not None else None
            # angles.append(cos_chi2) if cos_chi2 is not None else None
        return angles

    def _take_average_of_tuple(self, data):
        '''Takes an incoming tuple and computes the element-wise average of sets of lists, returning the result.'''
        result = ()
        for i in self.size:
            list_of_lists = []
            for j in self.size + 1:
                list_of_lists.append(data[i * (self.size + 1) + j])
            average = np.average(np.array(list_of_lists), axis=0)
            result += (average)
        return result