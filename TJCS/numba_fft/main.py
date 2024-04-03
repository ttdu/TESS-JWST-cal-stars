from numba_fft.data_gen import scramble_data
#from ls import lombscargle_fast
from numba_fft.fml import newfft

from numba import jit
import numpy as np


@jit(parallel=True, fastmath=True)
def false_alarms(datax, datay, iters=10000):
    ld = len(datax)
    maxes = np.array([newfft(*scramble_data(datax, datay, ld)) for i in range(iters)])
    # maxes = np.zeros(iters).astype(np.float64)
    # # Get an array of the maximum values of the periodogram
    # for i in range(iters):
    #     maxes[i] = np.nanmax(newfft(*scramble_data(datax, datay, ld)))
    # Get the 1,2,3 sigma values
    sigmas = np.quantile(maxes, [.682689, .954499, .997300])
    return sigmas