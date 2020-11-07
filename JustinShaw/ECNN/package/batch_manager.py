'''Formats Biopython's Structure object as keras-compatable dataset.'''

import numpy as np
from package.dwt import DWT as DWT
class BatchManager():

    def __init__(self, data, wavelet_size=4, verbose=False):
        '''Manages a moving window over the dataset based on the given wavelet size.
        
        Parameters
        ----------
        `data` - the biopython Structure object.
        `wavelet_size` - the input size to the dwt (must be power of two, defaults to 4)
        `verbose` - whether to display print statements.
        '''
        super().__init__()
        if verbose:
            print('Initializing Batch Manager...')
        self.data = data
        self.index = 0 # the left edge of the moving window
        self.num_angles = 44 # TODO: WTF is going on here????
        self.window_size = wavelet_size + wavelet_size * (wavelet_size + 1)
        
        # Check that the wavelet size is a power of two.
        if (wavelet_size & (wavelet_size - 1) == 0) and wavelet_size != 0:
            self.size = wavelet_size
        else:
            raise ValueError('wavelet_size must be a power of two!')

    def get_clean_dataset(self):
        '''Returns a cleaned dataset.'''
        X = []
        Y = []
        # Converts model-first data to residue-first data
        data_generator = self._next_window()
        for fine_data, coarse_data, output_data in data_generator:
            # Temporary storage for the DWT of the model-first data.
            fine_dwts = []
            coarse_dwts = []
            outputs = []

            # Convert data from model-first to angle-first grouping, then DWT the heck out of it.
            for grouped_fine_data in self._to_inverted_list(fine_data):
                fine_dwts.append(DWT(grouped_fine_data))
            for grouped_coarse_data in self._to_inverted_list(coarse_data):
                coarse_dwts.append(DWT(grouped_coarse_data))
            for grouped_output_data in self._to_inverted_list(output_data):
                outputs.append(np.array(grouped_output_data))
            
            # Loop over the number of phi/psi angles and add contents to X and Y.
            for i, fine_dwt in enumerate(fine_dwts):
                if len(X) == i: # create a new sublist to contain new inputs
                    X.append([])
                if len(Y) == i: # create a new sublist to contain new outputs
                    Y.append([])
                X[i].append(np.append(fine_dwt, coarse_dwts[i]))
                Y[i].append(np.array(outputs[i]))
        return np.array(X), np.array(Y)

    def _next_window(self):
        '''Generates the next `fine_data`, `coarse_data`, and `output_data` on the fly.
        
        ### Example:
        ```python
        >>> batch = BatchManager(data)
        >>> fine, coarse, output = batch.next()
        >>> fine # Returns a tuple of lists (wavelet_size)
        ([model_1], [model_2], [model_3], [model_4])
        >>> coarse # Returns a tuple of (wavelet_size) lists
        ([avg_1], [avg_2], [avg_3], [avg_4])
        >>> output # Returns a tuple of (wavelet_size + 1) lists
        ([model_1], [model_2], [model_3], [model_4], [model_5])
        >>> fine[0] # Returns a list of sin/cos phi/psi angles
        [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), sin(phi_2), cos(phi_2), ... , cos(psi_n)]
        >>> coarse[0] # Returns a list of the average of (wavelet_size + 1) sin/cos phi/psi angles
        [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), sin(phi_2), cos(phi_2), ... , cos(psi_n)]
        >>> output[0] # Returns a list of sin/cos phi/psi angles
        [sin(phi_1), cos(phi_1), sin(psi_1), cos(psi_1), sin(phi_2), cos(phi_2), ... , cos(psi_n)]
        ```
        '''
        # Yeild the next batch, but first check to make sure you're not out of data.
        while len(self.data) - self.index >= self.window_size + self.size + 1:
            try:
                # Define the range (slice) of data to take in the moving window.
                fine_slice = slice(self.index + self.window_size - self.size, self.index + self.window_size)
                coarse_slice = slice(self.index, self.index + self.window_size - self.size)
                output_slice = slice(self.index + self.window_size, self.index + self.window_size + self.size + 1)

                # Slice the data with the given slicing range.
                fine_data = self._slice_data(fine_slice)
                coarse_data = self._slice_data(coarse_slice, take_average=True)
                output_data = self._slice_data(output_slice)

                # Increment the moving window and return.
                self.index += 1
                yield fine_data, coarse_data, output_data
            except:
                raise StopIteration

    def _slice_data(self, data_slice, take_average=False):
        '''Returns a tuple of the list of angles in the given slice of the data.
        
        Parameters
        ----------
        `data_slice` - the slice of the data to parse.
        `take_average` - takes the average of the tuple if true (default false)
        '''
        result = ()
        models = self.data.child_list[data_slice]
        for model in models:
            angles = self._get_angles_from_model(model)
            result += tuple([angles])
        if take_average:
            result = self._take_average_of_tuple(result)
        return result

    def _get_angles_from_model(self, model, enhance_data=False):
        '''Returns a list of sin/cos phi/psi angles extracted from each model.
        
        Parameters
        ----------
        `model` - the model object to pull residues from.
        `enhance_data` - adds omega and chi2 angles to the list if true (default False)
        '''
        angles = []
        for residue in model.get_residues():
            self._add_angles_from_residue('phi', residue, angles)
            self._add_angles_from_residue('psi', residue, angles)
            if enhance_data:
                self._add_angles_from_residue('omega', residue, angles)
                self._add_angles_from_residue('chi2', residue, angles)
        while len(angles) < self.num_angles:
            angles.append(0)
        return angles
    
    def _add_angles_from_residue(self, angle_name, residue, angles_list):
        '''Appends the sin and cos of the given angle on the given residue to the angles list.
        
        Parameters
        ----------
        `angle_name` - the name of the angle to append.
        `residue` - the residue object search the dict for.
        `angles_list` - the array of angles to append.
        '''
        sin = residue.xtra.get(f'sin({angle_name})')
        cos = residue.xtra.get(f'cos({angle_name})')
        angles_list.append(sin) if sin is not None else None
        angles_list.append(cos) if cos is not None else None

    def _take_average_of_tuple(self, data):
        '''Takes an incoming tuple and computes the element-wise average of sets of lists, returning the result.
        
        Parameters
        ----------
        `data` - the large tuple of angle lists to condense.
        '''
        # Group the data into chunks of self.size (default 4)
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

    def _to_inverted_list(self, data):
        '''When given a tuple of lists, it returns an inverted list of lists.'''
        result = []
        for i in range(len(data)):
            for j in range(len(data[i])):
                if len(result) == j:
                    result.append([])
                result[j].append(data[i][j])
        return result