from numba import jit
import numpy as np

@jit 
def prod_ax0(array):
    """ Numba doesn't support the axis argument of np.prod
    """
    return np.array([np.prod(array[:, i]) for i in range(array.shape[1])])
    
@jit
def add_at(x, a, b):
    """ Replace the default numpy implementation
    """
    for i in np.arange(len(a)): x[a[i]] += b[i]

@jit
def extirpolate(x, y, N, M=4):
    """
    Extirpolate the values (x, y) onto an integer grid range(N),
    using lagrange polynomial weights on the M nearest points.

    Parameters
    ----------
    x : array-like
        array of abscissas
    y : array-like
        array of ordinates
    N : int
        number of integer bins to use. For best performance, N should be larger
        than the maximum of x
    M : int
        number of adjoining points on which to extirpolate.

    Returns
    -------
    yN : ndarray
        N extirpolated values associated with range(N)

    Examples
    --------
    >>> rng = np.random.default_rng(0)
    >>> x = 100 * rng.random(20)
    >>> y = np.sin(x)
    >>> y_hat = extirpolate(x, y)
    >>> x_hat = np.arange(len(y_hat))
    >>> f = lambda x: np.sin(x / 10)
    >>> np.allclose(np.sum(y * f(x)), np.sum(y_hat * f(x_hat)))
    True

    Notes
    -----
    This code is based on the C implementation of spread() presented in
    Numerical Recipes in C, Second Edition (Press et al. 1989; p.583).
    """
    #x, y = map(np.ravel, np.broadcast_arrays(x, y))


    #N = np.int64(np.max(x) + 0.5 * M + 1)


    # Now use legendre polynomial weights to populate the results array;
    # This is an efficient recursive implementation (See Press et al. 1989)
    result = np.zeros(N)

    # # first take care of the easy cases where x is an integer
    integers = x % 1 == 0
    add_at(result, x[integers].astype(np.int64), y[integers])
    
    x, y = x[~integers], y[~integers]
    

    # # For each remaining x, find the index describing the extirpolation range.
    # # i.e. ilo[i] < x[i] < ilo[i] + M with x[i] in the center,
    # # adjusted so that the limits are within the range 0...N
    ilo = np.clip((x - M // 2).astype(np.int64), 0, N - M)
    numerator = y * prod_ax0(x - ilo - np.arange(M)[:, np.newaxis])
    denominator = 6 #np.math.factorial(M - 1)
    

    for j in range(M):
        if j > 0:
            denominator *= j / (j - M)
        ind = ilo + (M - 1 - j)
        add_at(result, ind, numerator / (denominator * (x - ind)))
    return result