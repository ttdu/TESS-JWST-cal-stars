import numpy as np
import time
from astropy.timeseries import LombScargle

def simulate_pg(noise, size):
    # Generate a sample pg
    x = np.linspace(0, 28, size)
    y = np.random.normal(loc=0, scale=noise, size=size)
    
    # Get the standard deviation of y
    noise = np.full(size-2000, noise)
    
    # Delete a chunk from the middle
    x = np.concatenate((x[:10000], x[12000:]))
    y = np.concatenate((y[:10000], y[12000:]))
    
    # Get max/min periods
    max_f = 1/(2*(max(x)/size))
    min_f = 2/(max(x))
    df = max_f-min_f
    samples = int(5*max(x)*max_f)
    freqs = np.linspace(min_f, max_f, samples)
    
    # Create the periodogram
    LS = LombScargle(x, y, noise, normalization='psd')
    pwrs = LS.power(freqs)

    # Get the errors
    sigmas = np.array([.682689, .954499, .997300])
    alarms = LS.false_alarm_level(1-sigmas, method='bootstrap')
    return alarms

t0 = time.time()
simulate_pg(1, 20000)
print()