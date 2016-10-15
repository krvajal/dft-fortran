# this programs generate an interpolated 
# dataset for the chodrow potential
# as given in
# Burdick, G. A. (1963). 

# http://doi.org/10.1103/PhysRev.129.138

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
#load originla data file
data = np.loadtxt("copper_chodorow.txt",skiprows=1)
rmax = 2.4151
numpoints = 1000

x = np.linspace(0,rmax,numpoints);
pot =2 * 29*np.exp(-2.3151241717834*x**0.81266614122432+2.1984250222603e-02*x**4.2246376280056)-0.15595606773483*x-3.1350051440417e-03*x**2+5.1895222293006e-02*x**3 - 2.8027608685637e-02*x**4
plt.plot(data[:,0],data[:,1],'o',label="data")

f = interpolate.interp1d(data[:,0],data[:,1])



interpdata = np.zeros((numpoints,2))
interpdata[:,0] = x
interpdata[:,1] = f(x)

np.savetxt("linpot.txt",interpdata)

plt.plot(x,f(x),label='interpolation')
plt.plot(x,pot,label="param")
plt.legend()
plt.ioff()
plt.show()