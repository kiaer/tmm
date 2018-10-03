#advection diffusion model i 2 dimensioner
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


from sympy.abc import a, b

def velocity_field(psi): #insert stream function
    u = sympy.lambdify((a, b), psi.diff(b), 'numpy')
    v = sympy.lambdify((a, b), -psi.diff(a), 'numpy')
    return u, v


#parameters
diff    = 3.65 #diffusion constant [km^2 /y]
zmax    = 4 #depth of water column [m]
zGrid   = 20 #number of grid cells
#dz = zmax/(zGrid-1)
#z = np.arange(0 , zmax +dz, dz)
dz = zmax/(zGrid)
z = np.arange(0.5*dz , zmax, dz)
xmax    = 8000
xGrid   = 40
#dx = xmax/(xGrid-1)
#x = np.arange(0,xmax+dx,dx)
dx = xmax/(xGrid)
x = np.arange(0.5*dx,xmax,dx)
u = np.zeros((zGrid,xGrid))
w = np.zeros((zGrid,xGrid))
vmax = 4380 #[km/yr]

print('xvec_end=',x[-1])
print('zvec_end = ', z[-1])
print('shape x',np.shape(x))
print('shape z',np.shape(z))
print('dx= ',dx , dz)


def round_stream_function(A=vmax , H=zmax , L=xmax ): #L is rotation(lambda)
  return A*sympy.sin(a*(sympy.pi/L))*sympy.sin(b*(sympy.pi/H))

#extracting velocityfields from streamfunction and putting into arrays
psi = round_stream_function()
u_fun, w_fun = velocity_field(psi)

for i in range(0,zGrid):
    for j in range(0,xGrid):
        u[i,j] = u_fun(x[j],z[i])
        w[i,j] = w_fun(x[j],z[i])

print('maxw',np.amax(np.absolute(w)))
print('maxu',np.amax(np.absolute(u)))


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
C0[7,20] = 10  #initial conditions
#print("C0 =",C0)
C0f = C0.flatten()

#max time step for advection and diffusion
max_dtA = np.maximum(dx, dz)/np.maximum(np.amax(np.absolute(u)),np.amax(np.absolute(w)))
max_dtD = dz**2/(2*diff)

# choosing the minimum time step of the two
max_dt = np.minimum(max_dtA, max_dtD)

yrs = 5 #solve for number of years
yrs_out = np.arange(0,yrs,1) #output every time step

C = solve_ivp(lambda t, y: adflow2d(t, y, u, w, diff, dx, dz, xGrid, zGrid ), [0, yrs], C0f, t_eval=yrs_out, max_step=max_dt)
C.y = np.reshape(C.y,(zGrid,xGrid,len(C.t)))

concentration_sum = np.zeros(len(C.t))
for i in range(0,len(C.t)):
    concentration_sum[i] = np.sum(C.y[:,:,i])

plt.contourf(x,z, C0 ,cmap='jet')
plt.gca().invert_yaxis()
plt.xlabel('x')
plt.ylabel('depth')
plt.title('C0')
plt.colorbar()
print(plt.show())


plt.plot(C.t, concentration_sum)
plt.xlabel('t')
plt.ylabel('Sum of concentration')
plt.title('conservation test')
print(plt.show())

for i in range(1,10):
    plt.contourf(x,z, C.y[:,:,i] ,cmap='jet')
    plt.gca().invert_yaxis()
    plt.xlabel('x')
    plt.ylabel('depth')
    plt.title(['ode solution t=',i])
    plt.colorbar()
    print(plt.show())