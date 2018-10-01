import sympy
import matplotlib.pyplot as plt
import numpy as np

from sympy.abc import x, y

def round_stream_function(A=1, H = 10, L = 10): #L is rotation(lambda)
  return A*sympy.sin(sympy.pi * (x / L))*sympy.sin(sympy.pi * (y / H))#*sympy.sin(sympy.pi)

def round_stream_andy(A=1): #L is rotation(lambda)
  return A*sympy.sin(sympy.pi * sympy.cos(sympy.pi * x))*sympy.sin(sympy.pi*sympy.sin((sympy.pi*y) / 2))

def velocity_field(psi): #insert stream function
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


psi = round_stream_function()

u, w = velocity_field(psi)

ylim = (0, 10)
xlim = (0, 10)
fig, ax = plt.subplots(figsize=(4, 4))
plot_streamlines(ax, u, w, xlim, ylim)

c = plt.Circle((0, 0), radius=1, facecolor='none')
ax.add_patch(c)
print(sympy.sin(sympy.pi))
format_axes(ax)
print(u(10,10), w(10,10))
print(plt.show())
