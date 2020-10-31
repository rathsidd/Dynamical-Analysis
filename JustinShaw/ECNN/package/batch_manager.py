'''Generator that iterates and returns the next data for each batch.'''

import numpy as np

class BatchManager():

    def __init__(self, data, wavelet_size=4, verbose=False):
        '''Initialize the BatchManager class.
        
        This class manages a moving window over the given data based on the given wavelet size, and 
        '''
        super().__init__()
        if verbose:
            print('Initializing Batch Manager...')
        self.data = data
        self.index = 0
        self.num_angles = 44 # WTF
        self.size = wavelet_size # this should be a power of two
        self.window_size = wavelet_size + wavelet_size * (wavelet_size + 1)

    def next(self):
        '''Returns a tuple of the next `fine_data`, `coarse_data`, and `output_data`.'''
        fine_data, coarse_data = self._next_input()
        output_data = self._next_output()
        self.index += 1
        yield fine_data, coarse_data, output_data
    
    def _next_output(self):
        '''Returns a tuple of the next (default 5) data models.

        Returns
        -------
        The future data is returned first as a list of (default 5) raw elements.

        For example, this method might return:
        >>> ([model_1], [model_2], [model_3], [model_4], [model_5])
        where each model is a list that contains data, for example `model_5` might be:
        >>> [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), sin(phi_2), cos(phi_2), ... , cos(psi_n)]
        '''
        try:
            if len(self.data) - self.index >= self.window_size + self.size + 1:
                y_train = self._get_output_data()
                return y_train
        except:
            raise StopIteration
    
    def _next_input(self):
        '''Returns a tuple of fine and coarse data.
        
        Returns
        -------
        The fine data is returned first as a list of (default 4) raw elements.
        The coarse data is returned second as a list of (default 4) averages of (default 5) elements.
        
        For example, this method might return:
        >>> ([model_1], [model_2], [model_3], [model_4]), ([avg_5], [avg_6], [avg_7], [avg_8])
        where each model is a list that contains data, for example `model_4` might be:
        >>> [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), sin(phi_2), cos(phi_2), ... , cos(psi_n)]
        '''
        try:
            if len(self.data) - self.index >= self.window_size:
                fine_data = self._get_fine_data()
                coarse_data = self._get_coarse_data()
                return fine_data, coarse_data
        except:
            raise StopIteration

    def _get_output_data(self):
        '''Returns a tuple of length `size + 1` with elements that are themselves lists of sin/cos phi/psi angles.'''
        result = ()
        models = self.data.child_list[self.index + self.window_size : self.index + self.window_size + self.size + 1]
        for model in models:
            angles = self._get_angles_from_model(model)
            result += tuple([angles])
        return result

    def _get_fine_data(self):
        '''Returns a tuple of length `size` with elements that are themselves lists of sin/cos phi/psi angles.'''
        result = ()
        models = self.data.child_list[self.index + self.window_size - self.size : self.index + self.window_size]
        for model in models:
            angles = self._get_angles_from_model(model)
            result += tuple([angles])
        return result

    def _get_coarse_data(self):
        '''Returns a tuple of length `size` with elements that are themselves lists of the average of sin/cos phi/psi angles.'''
        result = ()
        models = self.data.child_list[self.index : self.index + self.window_size - self.size]
        for model in models:
            angles = self._get_angles_from_model(model)
            result += tuple([angles])
        return self._take_average_of_tuple(result)

    def _get_angles_from_model(self, model):
        '''Returns a list of sin/cos phi/psi angles extracted from each model.'''
        angles = []
        for residue in model.get_residues():
            sin_phi = residue.xtra.get('sin(phi)')
            cos_phi = residue.xtra.get('cos(phi)')
            sin_psi = residue.xtra.get('sin(psi)')
            cos_psi = residue.xtra.get('cos(psi)')
            # sin_omega = residue.xtra.get('sin(omega)')
            # cos_omega = residue.xtra.get('cos(omega)')
            # sin_chi2 = residue.xtra.get('sin(chi2)')
            # cos_chi2 = residue.xtra.get('cos(chi2)')
            angles.append(sin_phi) if sin_phi is not None else None
            angles.append(cos_phi) if cos_phi is not None else None
            angles.append(sin_psi) if sin_psi is not None else None
            angles.append(cos_psi) if cos_psi is not None else None
            # angles.append(sin_omega) if sin_omega is not None else None
            # angles.append(cos_omega) if cos_omega is not None else None
            # angles.append(sin_chi2) if sin_chi2 is not None else None
            # angles.append(cos_chi2) if cos_chi2 is not None else None
        while len(angles) < self.num_angles:
            angles.append(0)
        return angles

    def _take_average_of_tuple(self, data):
        '''Takes an incoming tuple and computes the element-wise average of sets of lists, returning the result.'''
        result = ()
        for i in range(self.size):
            list_of_lists = np.zeros(shape=(self.size + 1, len(data[0])))
            for j in range(self.size + 1):
                sub_list = data[i * (self.size + 1) + j]
                while len(sub_list) < self.num_angles:
                    sub_list.append(0)
                list_of_lists[j] = sub_list
            average = np.average(list_of_lists, axis=0)
            result += tuple([average])
        return result