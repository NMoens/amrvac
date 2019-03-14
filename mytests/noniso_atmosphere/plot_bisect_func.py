import numpy as np
import matplotlib.pyplot as plt


c0 = -112909.69294058459
c1 = 4.4770603528486008

def poly(x):
    return x**4.0 + x*c1 - c0


xx = np.linspace(-abs(c0/c1),abs(c0/c1),1000)

plt.plot(xx,poly(xx))
plt.show()
