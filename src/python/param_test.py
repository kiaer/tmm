#parameters for advection diffusion model

import numpy as np

diff = 4.3 #diffusion constant [m^2 /d]
u = 1 # sinking velocity [m/d]
zmax = 200 #depth of water column [m]
dz = 0.5 #grid size
nGrid = int(zmax/dz) #number of grid cells
z = np.mgrid[0.5*dz:zmax:dz] # depth vector located in the middle of each grid cell

C0 = z*0  #initial cond
