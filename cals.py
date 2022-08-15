import numpy as np


g=9.80665
k=0.0076
deep=17.56

freq = np.sqrt(g * k * np.tanh(k * deep))

print(freq)