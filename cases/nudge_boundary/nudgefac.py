import matplotlib.pyplot as pl
import numpy as np

pl.close('all')

# Numerical grid
xsize = 51200
dx    = 50
x     = np.arange(0.5*dx, xsize, dx)

# Boundary settings
dxb   = 300     # Width of boundary nudging
oxb   = 1025    # Offset from lateral boundary

# Center of nudging area
bc1 = oxb
bc2 = xsize-oxb

f1 = np.exp(-((x-bc1)**2/(2*dxb**2)))
f2 = np.exp(-((x-bc2)**2/(2*dxb**2)))

pl.figure()
pl.plot(x, f1)
pl.plot(x, f2)
