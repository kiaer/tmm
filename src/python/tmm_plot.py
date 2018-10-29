from typing import Any, Union

import matplotlib.pyplot as plt
# import param_test #parameters
import numpy as np
from numpy.linalg import matrix_power
from scipy.integrate import solve_ivp
import sympy
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix
import scipy
import scipy.io as io
from sympy.abc import a, b

diff    = 1 #diffusion constant [km^2 /y]
zmax    = 50 #depth of water column [m]
zGrid   = 50 #number of grid cells
dz = zmax/(zGrid-1)
z = np.arange(0 , zmax +dz, dz)
xmax    = 100
xGrid   = 50
dx = xmax/(xGrid-1)
x = np.arange(0,xmax+dx,dx)
u = np.zeros((zGrid,xGrid))
w = np.zeros((zGrid,xGrid))

A = scipy.sparse.load_npz('../../bin/A.npz')
#io.savemat('A.mat', {'A': A})

print(np.shape(A))
print(A)

plt.contourf(A.todense(), cmap='jet')
plt.gca().invert_yaxis()
plt.colorbar()
print(plt.show())

plt.spy(A, marker=None, markersize=None)
plt.gca().invert_yaxis()
plt.show()



Apow = np.linalg.matrix_power(A.todense(), 500)
plt.contourf(Apow, cmap='jet')
plt.gca().invert_yaxis()
plt.colorbar()
print(plt.show())
Apow1 = A ** 500

C0 = np.zeros((zGrid,xGrid))  #initial condition matrix
C0[10,10] = 1  #initial conditions
#print("C0 =",C0)
C0f = C0.flatten()

Ca10f = Apow.dot(C0f)
Ca10 = np.reshape(Ca10f,(zGrid,xGrid))

plt.contourf(x,z, Ca10 ,cmap='jet')
plt.gca().invert_yaxis()
plt.xlabel('x')
plt.ylabel('depth')
plt.title('ode solution t=10')
plt.colorbar()
print(plt.show())




