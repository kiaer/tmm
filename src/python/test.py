import numpy as np
from sympy import sin, pi

a = np.zeros(8)
a[0:2] = 1
a[2:4] = 2
a[4:6] = 3
a[6:8] = 4


b = a[0:len(a):2]
print(b)

print(sin(pi))
