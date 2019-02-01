import numpy as np
import matplotlib.pyplot as plt


c0 = 5.0179671072039653
c1 = 1.0000352601828886

def poly(x):
    return x**4.0 + x*c1 - c0


xx = np.linspace(0,abs(c0/c1),100)

print poly(xx)

plt.plot(xx,poly(xx))
plt.show()
