"""Calculate the electric field due to a dipole.
"""
import numpy as np
import matplotlib.pyplot as plt

# Don't want to rely on anyone having SciPy installed as we only want it for
# the definition of the permittivity of free space, if it's not installed then
# we just define it here.
try:
	import scipy.constants.epsilon_0 as epsilon_0
except ImportError:
	epsilon_0 = 8.854e-12

import util

# Side dimension [m] of the dipole and magnitude [C] of the charges.
d = 1.0
qmag = 1.0

# Construct a 2D grid.
xe, ye, xc, yc = util.cartesian_2d_grid(-5.15, 5.15, -10.15, 10.15, 0.04, 0.04)

# Figure 20cm x 10cm (the dimensions in Matplotlib are in inches, 1"=2.54 cm).
fig = plt.figure(figsize=(20/2.54,10/2.54))

# Calculate the potential, then take the grad to get the electric field.  Note
# that if the grid resolution is changed above then it will need to be changed
# here too.
V = util.Vmono(-qmag, xc, yc-0.5*d) + util.Vmono(qmag, xc, yc+0.5*d)
Ex = -np.gradient(V, 0.04, axis=1)
Ey = -np.gradient(V, 0.04, axis=0)

# Make the plot.
hc = plt.contourf(xc, yc, V, np.linspace(-1e11,1e11,100), cmap='RdBu_r', extend='both')
plt.contour(xc, yc, V, [-1e10,-1e9,-1e8,-1e7,1e7,1e8,1e9,1e10], colors='k')
plt.axis('equal')
plt.streamplot(xc, yc, Ex, Ey, color='#555555',linewidth=0.5)
plt.xlim([-6,6])
plt.title('Electric dipole')
hcb = plt.colorbar(hc)
hcb.set_label('Potential [V]')

plt.show()
