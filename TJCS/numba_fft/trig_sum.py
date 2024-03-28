from numba import jit
import numpy as np

from numba_fft.extirpolate import extirpolate

@jit
def my_bitceil(N):
    return np.int64(2 ** np.ceil(np.log2(N)))

@jit
def trig_sum(t, h, df, N, f0=0, freq_factor=1, oversampling=5, Mfft=4):
    """Compute (approximate) trigonometric sums for a number of frequencies.

    This routine computes weighted sine and cosine sums::

        S_j = sum_i { h_i * sin(2 pi * f_j * t_i) }
        C_j = sum_i { h_i * cos(2 pi * f_j * t_i) }

    Where f_j = freq_factor * (f0 + j * df) for the values j in 1 ... N.
    The sums can be computed either by a brute force O[N^2] method, or
    by an FFT-based O[Nlog(N)] method.

    Parameters
    ----------
    t : array-like
        array of input times
    h : array-like
        array weights for the sum
    df : float
        frequency spacing
    N : int
        number of frequency bins to return
    f0 : float, optional
        The low frequency to use
    freq_factor : float, optional
        Factor which multiplies the frequency
    use_fft : bool
        if True, use the approximate FFT algorithm to compute the result.
        This uses the FFT with Press & Rybicki's Lagrangian extirpolation.
    oversampling : int (default = 5)
        oversampling freq_factor for the approximation; roughly the number of
        time samples across the highest-frequency sinusoid. This parameter
        contains the trade-off between accuracy and speed. Not referenced
        if use_fft is False.
    Mfft : int
        The number of adjacent points to use in the FFT approximation.
        Not referenced if use_fft is False.

    Returns
    -------
    S, C : ndarray
        summation arrays for frequencies f = df * np.arange(1, N + 1)
    """
    df *= freq_factor
    f0 *= freq_factor

    if df <= 0:
        raise ValueError("df must be positive")
    #t, h = map(np.ravel, np.broadcast_arrays(t, h))

    #if use_fft:
    Mfft = int(Mfft)
    if Mfft <= 0:
        raise ValueError("Mfft must be positive")

    # required size of fft is the power of 2 above the oversampling rate
    Nfft = my_bitceil(N * oversampling)
    t0 = np.min(t)

    if f0 > 0:
        # numba is unhappy without this
        j = h.astype(np.complex128)
        j = j * np.exp(2j*np.pi*f0*(t-t0))
    #return j
    #print("prod")
    #print(np.mod((t-t0)*df*Nfft, Nfft))
    tnorm = np.mod((t-t0)*df*Nfft, Nfft)
    grid = extirpolate(tnorm, h, Nfft, Mfft)
    
    fftgrid = np.fft.ifft(grid)[:N]
    
    
    # # t0 should always be zero
    # # if t0 != 0:
    # #     f = f0 + df * np.arange(N)
    # #     fftgrid *= np.exp(2j * np.pi * t0 * f)

    C = Nfft * fftgrid.real
    S = Nfft * fftgrid.imag

    return S, C