import sympy
import matplotlib.pyplot as plt
import numpy as np

from sympy.abc import x, y


def cylinder_stream_function(U=1, R=5):
    r = sympy.sqrt(x**2 + y**2)
    theta = sympy.atan2(y, x)
    #print(theta)
    return U * (r - R**2 / r) * sympy.sin(theta)
    #return U * sympy.sin(theta)

def corner_stream_function(n=-3, A=1):
  r = sympy.sqrt(x**2 + y**2)
  theta = sympy.atan2(y, x)
  return A * r**n * sympy.sin(n * theta)

# Starter her på egen stream
# ??
def round_stream_function(A=2, H = 5, L = 3, n = -2):
  r = sympy.sqrt(x**2 + y**2)
  theta = sympy.atan2(y, x)
  return r ** n * sympy.cos((x / L) + (y / H))
  #return r ** n * sympy.sin(theta) 


def velocity_field(psi):
    u = sympy.lambdify((x, y), psi.diff(y), 'numpy')
    v = sympy.lambdify((x, y), -psi.diff(x), 'numpy')
    return u, v

def plot_streamlines(ax, u, v, xlim=(-1, 1), ylim=(-1, 1)):
    x0, x1 = xlim
    y0, y1 = ylim
    Y, X =  np.ogrid[y0:y1:100j, x0:x1:100j]
    ax.streamplot(X, Y, u(X, Y), v(X, Y), color='cornflowerblue')

def format_axes(ax):
    ax.set_aspect('equal')
    ax.figure.subplots_adjust(bottom=0, top=1, left=0, right=1)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    #for spine in ax.spines.items():
    #    spine.set_visible(False)


#psi = cylinder_stream_function()
#psi = corner_stream_function()
psi = round_stream_function()
print(psi)

u, v = velocity_field(psi)

print(u(0,0))
print(v(1,1))

xlim = ylim = (-3, 3)
fig, ax = plt.subplots(figsize=(4, 4))
plot_streamlines(ax, u, v, xlim, ylim)

c = plt.Circle((0, 0), radius=1, facecolor='none')
ax.add_patch(c)

format_axes(ax)

print(plt.show())
