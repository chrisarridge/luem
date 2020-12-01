"""Calculate the magnetic field from a simple solar source surface model.
"""
import numpy as np
import matplotlib.pyplot as plt

import util

# Radius of the source surface and strength of the dipole.
R = 10.0
k = 1.0

# Construct a 2D grid.
xe, ze, xc, zc = util.cartesian_2d_grid(-15, 15, -15, 15, 0.25, 0.25)

# Figure 20cm x 10cm (the dimensions in Matplotlib are in inches, 1"=2.54 cm).
fig = plt.figure(figsize=(20/2.54,10/2.54))

# Compute spherical coordinates and calculate the field.
r = np.sqrt(xc**2 + zc**2)
ph = np.arctan2(0.0, xc)
th = np.arccos(zc/r)
Br = (k/r**3)*(2 + (r/R)**3)*np.cos(th)
Bth = (k/r**3)*(1 - (r/R)**3)*np.sin(th)
Bth[r>R] = 0.0
Bx = (Br*np.sin(th) + Bth*np.cos(th))*np.cos(ph)
Bz = Br*np.cos(th) - Bth*np.sin(th)

plt.streamplot(xc, zc, Bx, Bz)
plt.plot(R*np.cos(np.linspace(0,2*np.pi,64)),R*np.sin(np.linspace(0,2*np.pi,64)),'-r')
plt.axis('equal')
#
#
# # Make the plot.
# hc = plt.contourf(xc, yc, V, np.linspace(-1e11,1e11,100), cmap='RdBu_r', extend='both')
# plt.contour(xc, yc, V, [-1e10,-1e9,-1e8,-1e7,1e7,1e8,1e9,1e10], colors='k')
# plt.axis('equal')
# plt.streamplot(xc, yc, Ex, Ey, color='#555555',linewidth=0.5)
# plt.xlim([-6,6])
# plt.title('Electric dipole')
# hcb = plt.colorbar(hc)
# hcb.set_label('Potential [V]')

plt.show()
