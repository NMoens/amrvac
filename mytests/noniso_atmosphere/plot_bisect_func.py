import numpy as np
import matplotlib.pyplot as plt


c0 = 7.0e73
c1 = 1.4e55

def poly(x):
    return x**4.0 + x*c1 - c0


xx = np.linspace(-abs(c0/c1),abs(c0/c1),100)

plt.plot(xx,poly(xx))
plt.show()
