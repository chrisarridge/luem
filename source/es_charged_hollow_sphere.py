"""This code calculates the electric field due to a charged hollow sphere.
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


# The total charge and radius of the sphere.
R = 1.0
Q = 1.0

# Setup a grid of surface area patches on the sphere with an integer number
# of patches in polar angle (theta) and azimuthal angle (phi).  The _th and
# _ph positions are the centres of each patch.
nth = 100
dth = np.pi/nth
_th = np.linspace(0, np.pi, (nth+1))[0:-1] + dth/2.0

nph = 200
dph = 2*np.pi/nph
_ph = np.linspace(0, 2*np.pi, (nph+1))[0:-1] + dph/2.0

th, ph = np.meshgrid(_th, _ph)
th = th.flatten()
ph = ph.flatten()

# Calculate the position of each patch in Cartesian coordinates.
xp = R*np.cos(ph)*np.sin(th)
yp = R*np.sin(ph)*np.sin(th)
zp = R*np.cos(th)

# Calculate the area of each patch.
dS = R*R*np.sin(th)*dth*dph

# Calculate the charge on each patch from the total charge, the area of a sphere
# to calculate the surface charge density, and the area of the patch.
dq = dS*Q/(4*np.pi*R*R)

# Calculate the field on the x axis as a function of x.
x = np.linspace(0, 10, 101)
y = np.zeros_like(x)
z = np.zeros_like(x)

# Clear the electric field to start the summation.
Ex = np.zeros(len(x))
Ey = np.zeros(len(x))
Ez = np.zeros(len(x))

util.do_patches(x, y, z, xp, yp, zp, dq, Ex, Ey, Ez)

# Add physical constants.
Ex *= 1.0/(4.0*np.pi*epsilon_0)
Ey *= 1.0/(4.0*np.pi*epsilon_0)
Ez *= 1.0/(4.0*np.pi*epsilon_0)

# Plot the field strength in each component as a function of x.
plt.plot(x, Ex, label=r'$E_x$')
plt.plot(x, Ey, label=r'$E_y$')
plt.plot(x, Ez, label=r'$E_z$')
plt.legend()
plt.xlabel('x [m]')
plt.ylabel('Electric field [V/m]')
plt.show(block=True)
