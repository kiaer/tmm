#advection diffusion model i 2 dimensioner
from typing import Any, Union

import matplotlib.pyplot as plt
# import param_test #parameters
import numpy as np
from numpy.linalg import matrix_power
from scipy.integrate import solve_ivp
import sympy


from sympy.abc import a, b

def velocity_field(psi): #insert stream function
    u = sympy.lambdify((a, b), psi.diff(b), 'numpy')
    v = sympy.lambdify((a, b), -psi.diff(a), 'numpy')
    return u, v


#parameters
diff    = 0.1 #diffusion constant [m^2 /d]
zmax    = 10 #depth of water column [m]
zGrid   = 30 #number of grid cells
dz = zmax/(zGrid-1)
z = np.mgrid[0.5*dz : zmax+dz : dz]
xmax    = 10
xGrid   = 30
dx = xmax/(xGrid-1)
x = np.mgrid[0.5*dx:xmax+dx:dx]
u = np.zeros((zGrid,xGrid))
w = np.zeros((zGrid,xGrid))

print('xvec_end=',x[-1])
print('zvec_end = ', z[-1])
print('shape x',np.shape(x))


def round_stream_function(A=1 , H=zmax , L=xmax ): #L is rotation(lambda)
  return A*sympy.sin(a*(sympy.pi/L))*sympy.sin(b*(sympy.pi/H))

#extracting velocityfields from streamfunction and putting into arrays
psi = round_stream_function()
u_fun, w_fun = velocity_field(psi)

for i in range(0,zGrid):
    for j in range(0,xGrid):
        u[i,j] = u_fun(x[j],z[i])
        w[i,j] = w_fun(x[j],z[i])



# the adv-diff model
def adflow2d(t, c, u, w, diff, dx, dz, xGrid, zGrid):
    #advective flux, xdir
    c = np.reshape(c,(zGrid,xGrid))
    #print(np.shape(c))
    Jax = np.zeros((zGrid,xGrid+1))
    Jaz = np.zeros((zGrid + 1, xGrid))
    Jdx = np.zeros((zGrid,xGrid+1))
    Jdz = np.zeros((zGrid + 1, xGrid))
    for k in range(1,zGrid):
        for l in range(1, xGrid):
            if u[k,l] >= 0:
                Jax[k, l] = u[k,l-1] * c[k,l-1]
            else:
                Jax[k, l] = u[k,l]*c[k,l]

            # advective flux, zdir
            if w[k,l] >= 0:
                Jaz[k,l] = w[k-1,l] * c[k-1 , l]
            else:
                Jaz[k, l] = w[k, l] * c[k, l]

            Jaz[0, l] = 0
            Jaz[-1, l] = 0

        Jax[k, 0] = 0
        Jax[k, -1] = 0
    for m in range(0,zGrid):
        # diff flux x-dir
        Jdx[m, 1:xGrid] = -diff * ((c[m, 1:xGrid] - c[m, 0:xGrid - 1]) / dx)
        # diff flux z dir
    for n in range(1,xGrid):
        Jdz[1:zGrid, n] = -diff * ((c[1:zGrid, n] - c[0:zGrid-1 , n]) / dz)
    Jdx[:,0] = 0
    Jdx[:,-1] = 0
    Jdz[0,:] = 0
    Jdz[-1,:] = 0

    #adv+diff
    Jx = Jax #+ Jdx
    Jz = Jaz + Jdz
    dcdt = (Jx[:,0:xGrid]-Jx[:,1:xGrid+1])/dx + (Jz[0:zGrid,:]-Jz[1:zGrid+1,:])/dz
    dcdtf = dcdt.flatten()
    return dcdtf


C0 = np.zeros((zGrid,xGrid))  #initial condition matrix
C0[5:7,8:10] = 1  #initial conditions
#print("C0 =",C0)
C0f = C0.flatten()

#max time step for advection and diffusion
max_dtA = np.maximum(dx, dz)/np.maximum(np.amax(np.absolute(u)),np.amax(np.absolute(w)))
max_dtD = dz**2/(2*diff)

# choosing the minimum time step of the two
max_dt = np.minimum(max_dtA, max_dtD)


C = solve_ivp(lambda t, y: adflow2d(t, y, u, w, diff, dx, dz, xGrid, zGrid ), [0, 50], C0f, max_step=max_dt) #kan tilføje , max_step =
C.y = np.reshape(C.y,(zGrid,xGrid,len(C.t)))

print('ode done')
#building transport matrix, A
A = np.zeros((xGrid*zGrid,xGrid*zGrid))
A_ind = 0
Ca0 = np.zeros((zGrid,xGrid))  #initial condition matrix

for i in range(0,zGrid):
    for j in range(0,xGrid):
        Ca0[i,j] = 1 #insert tracer
        Ca0f = Ca0.flatten()
        Ca = solve_ivp(lambda t, y: adflow2d(t, y, u, w, diff, dx, dz, xGrid, zGrid ) , [0, 1], Ca0f, max_step=max_dt)
        Ca.y = np.reshape(Ca.y,(zGrid,xGrid,len(Ca.t)))
        A_loc = Ca.y[:,:,-1]
        A[:,A_ind] = A_loc.flatten()
        Ca0[i,j] = 0 #remove tracer
        A_ind = A_ind +1


from numpy.linalg import matrix_power

Ca10 = np.zeros((zGrid,xGrid))
An = matrix_power(A,50)

print('shape.C0f',np.shape(C0f))
print('shape.C0f Tr',np.shape(C0f.transpose()))
print('A^10',np.shape(An))

Ca10f = An.dot(C0f)

Ca10 = np.reshape(Ca10f,(zGrid,xGrid))






print('shape C.t',len(C.t))
print('shape C.y',np.shape(C.y[:,:,-1]))
print('shape Ca10',np.shape(Ca10))


concentration_sum = np.zeros(len(C.t))
for i in range(0,len(C.t)):
    concentration_sum[i] = np.sum(C.y[:,:,i])

plt.plot(C.t, concentration_sum)
plt.xlabel('t')
plt.ylabel('Sum of concentration')
plt.title('conservation test')
print(plt.show())

#print(np.shape(C.y))
#print("C = ",C.y[:,:,-1])
#for i in range(0,len(C.t),5):
 #   plt.contourf(x,z, C.y[:,:,i] ,cmap='jet')
  #  plt.gca().invert_yaxis()
   # plt.xlabel('x')
#    plt.ylabel('depth')
 #   plt.title(['ode solution t=',str(C.t[i])])
  #  plt.colorbar()
   # print(plt.show())

plt.contourf(x,z, C.y[:,:,-1] ,cmap='jet')
plt.gca().invert_yaxis()
plt.xlabel('x')
plt.ylabel('depth')
plt.title('ode solution t=10')
plt.colorbar()
print(plt.show())



plt.contourf(x,z, Ca10 ,cmap='jet')
plt.gca().invert_yaxis()
plt.xlabel('x')
plt.ylabel('depth')
plt.title('matrix solution A t= 10')
plt.colorbar()
print(plt.show())


