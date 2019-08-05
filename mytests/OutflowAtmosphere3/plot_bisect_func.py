import numpy as np
import matplotlib.pyplot as plt


c0 = 3124.0
c1 = 1.000

def poly(x):
    return x**4.0 + x*c1 - c0

up_x = max(abs(c0/c1), abs(c0)**(1./4.))
xx = np.linspace(-up_x,up_x,10000000)

plt.plot(xx,poly(xx))
plt.plot(xx,poly(xx)*0.0)
plt.show()
