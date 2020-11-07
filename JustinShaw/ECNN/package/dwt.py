import pywt
import numpy as np

def DWT(self, data):
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
        (cA, cD) = pywt.dwt(data, 'haar') # TODO Is this the right transform?
        return np.append(cA, cD)
    else:
        raise ValueError('Data should contain a power of two number of elements')