from numba import jit
import numpy as np

@jit()
def scramble_data(datax, datay, ld):
    # Generate random integers for indices of data
    s = np.random.randint(0, ld, ld)
    return datax, datay[s], ld
    
    # for _ in range(n_bootstrap):
    #     s = rng.integers(0, len(y), len(y))  
    #     ls_boot = LombScargle(
    #         t, y[s], dy if dy is None else dy[s], normalization=normalization
    #     )
    #     freq, power = ls_boot.autopower(maximum_frequency=fmax)
    #     power_max.append(power.max())
    # x = np.linspace(0, 28, size)
    # y = np.random.normal(loc=0, scale=noise, size=size)
    # y += 1.7*np.sin(np.pi*x)
    
    # # Delete a chunk from the middle
    # x = np.concatenate((x[:10000], x[12000:]))
    # y = np.concatenate((y[:10000], y[12000:]))
    
    # # Get the standard deviation of y
    # noise = np.full(size-2000, noise).astype(np.float64)
    
    # # Get max/min periods
    # max_f = 1/(2*(max(x)/size))
    # min_f = 2/(max(x))
    # df = max_f-min_f
    # Nsamples = int(max(x)*max_f)

    # return x, y, noise, min_f, df, Nsamples