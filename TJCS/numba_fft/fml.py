from numba import jit
import numpy as np
from numba_fft.trig_sum import trig_sum

@jit(fastmath=True)
def newfft(x, y, ld):
    """Inputs:
    x: value of "x", in this case time
    y: value of "y", in this case flux
    ld: length of data
    """
    
    # Nyquist frequencies
    max_f = 1/(2*(max(x)/ld))
    min_f = 2/(max(x))
    df = max_f-min_f
    Nsamples = int(max(x)*max_f)

    # w = np.power(noise, -2)
    # w /= np.sum(w)
    # Since we're setting noise to 1, we can just set this value
    w = np.full(len(y), 1/len(y)).astype(np.float64)
    
    y -= np.dot(w, y)


    Sh, Ch = trig_sum(t=x, h=w*y, f0=min_f, df=df, N=Nsamples)
    S2, C2 = trig_sum(t=x, h=w, freq_factor=2, f0=min_f, df=df, N=Nsamples)

    S, C = trig_sum(x, w, f0=min_f, df=df, N=Nsamples)
    tan_2omega_tau = (S2 - 2 * S * C) / (C2 - (C * C - S * S))
    

    S2w = tan_2omega_tau / np.sqrt(1 + tan_2omega_tau * tan_2omega_tau)
    C2w = 1 / np.sqrt(1 + tan_2omega_tau * tan_2omega_tau)
    Cw = np.sqrt(0.5) * np.sqrt(1 + C2w)
    Sw = np.sqrt(0.5) * np.sign(S2w) * np.sqrt(1 - C2w)

    #YY = np.dot(w, y**2)
    YC = Ch * Cw + Sh * Sw
    YS = Sh * Cw - Ch * Sw
    CC = 0.5 * (1 + C2 * C2w + S2 * S2w)
    SS = 0.5 * (1 - C2 * C2w - S2 * S2w)

    CC -= (C * Cw + S * Sw) ** 2
    SS -= (S * Cw - C * Sw) ** 2

    # return the amplitude
    amp = np.sqrt(2*(YC * YC / CC + YS * YS / SS))

    # consequences of fast_math...
    newamp = amp[np.argwhere(amp!=np.inf).flatten()]
    
    return np.nanmax(newamp)