from numba_fft.fml import newfft
from numba_fft.fml import trig_sum
import numpy as np

x = np.linspace(1, 10000, 10000)
y = np.random.randint(1, 100, 10000)*6.28
l = len(y)

for _ in range(10000):
    newfft(x, y, l)

# for _ in range(200000):
#     trig_sum(x, y, 0.01, 20000, 0.01)