#advection diffusion model i 1 dimension

import param_test #parameters
import numpy as np
import matplotlib.pyplot as plt

from numpy.linalg import matrix_power
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

def adflow1d2(t, c, u, diff, dz, nGrid):
    #advective flux
    Ja = np.zeros(nGrid+1)
    Ja[1:nGrid] = u * c[0:nGrid-1]
    Ja[0] = 0
    Ja[-1] = 0
    #diffusive flux
    Jd = np.zeros(nGrid+1)
    Jd[1:nGrid] = -diff * ((c[1:nGrid]-c[0:nGrid-1])/dz)
    Jd[0] = 0
    Jd[-1] = 0
    #adv+diff
    J = Ja + Jd
    dcdt = (J[0:nGrid]-J[1:nGrid+1])/dz
    return dcdt

time_vec = np.linspace(0 , 1 , 30)
C0 = param_test.C0 #initial conditions
s = param_test.nGrid
A = np.zeros((s,s))
for i in range(len(C0)):  #len(C0)
    C0[i] = 1 #insert tracer
    C = solve_ivp(lambda t, y: adflow1d2(t, y, param_test.u, param_test.diff, param_test.dz, param_test.nGrid) , [0, 1], C0, t_eval=[1])
    A[:,i] = C.y[:,-1]
    C0[i] = 0 #remove tracer



#print(np.ndim(C))

#print("C=", C.shape)
#print("time = ",time_vec.shape)
np.savetxt('A_mat.txt',A,delimiter=',',newline=';')

time_vec2 = np.linspace(0 , 12 , 160)
C0[int(len(C0) / 2)] = 1 #insert tracer
#print(C0)
#C = odeint(adflow1d,C0,time_vec2,args = (param_test.u,param_test.diff,param_test.dz, param_test.nGrid))

C = solve_ivp(lambda t, y: adflow1d2(t, y, param_test.u, param_test.diff, param_test.dz, param_test.nGrid) , [0, 30], C0)

#print(C.y[:,-1])
#print(C.t.shape)


#plt.contourf(time_vec,param_test.z, C.T ,cmap='jet')
#plt.gca().invert_yaxis()
#plt.xlabel('days')
#plt.ylabel('depth')
#plt.title('')
#print(plt.show())
plt.figure(1)
plt.subplot(1,3,1)
plt.gca().invert_yaxis()
plt.plot(C.y[:,-1],param_test.z)
plt.title('ode solution')

C2 = matrix_power(A,30).dot(C0)

plt.subplot(1,3,2)
plt.gca().invert_yaxis()
plt.plot(C2,param_test.z)
plt.title('matrix solution')

delta = C.y[:,-1]-C2

plt.subplot(1,3,3)
plt.gca().invert_yaxis()
plt.plot(delta,param_test.z)
plt.title('ODE - matrix')
print(plt.show())

np.set_printoptions(precision=3)
np.set_printoptions(threshold=np.nan)

print(A)

