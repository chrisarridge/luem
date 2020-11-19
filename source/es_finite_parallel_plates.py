"""This code calculates the electric field due to a pair of finite parallel plates.
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



def charged_unit_parallelplates(x, y, z, Qtop, Qbottom, d, n=10):
	"""Calculates the electric field of a unit parallel plate.

	The plates are split into n^2 small rectangles (patches) and the
	total charge on each plate is divided up into these small rectangles.  Then
	the Coulomb law is used to sum up the field contributions from all the
	patches.  The plates are assumed to lie in the x-y plane a distance z=+d/2
	and z=-d/2.

	Args:

		:param x: X coordinates for the points to calculate the field at [m].
		:param y: X coordinates for the points to calculate the field at [m].
		:param z: Z coordinates for the points to calculate the field at [m].
		:param Qtop: Total charge on the top (+d/2) plate [C].
		:param Qbottom: Total charge on the bottom (-d/2) plate [C].
		:param n: Optional, number of cells to split in each dimension (default=10).
		:return: Ex [V/m], Ey [V/m], Ez [V/m]
	"""

	# Setup a grid of points over the face of a unit plane (x=-0.5->x=0.5,
	# y=-0.5->y=0.5); and the x and y axes are split into n chunks between
	# -0.5 and 0.5, producing n^2 patches over the face.
	u, v = np.meshgrid(np.linspace(-0.5,0.5,n), np.linspace(-0.5,0.5,n))
	dp = u[0,1]-u[0,0]
	dS = dp*dp		# area of the patch - these are all uniform.
	u = u.flatten()
	v = v.flatten()

	# Clear the electric field to start the summation.
	Ex = np.zeros(len(x))
	Ey = np.zeros(len(x))
	Ez = np.zeros(len(x))

	# The surface charge density is just the total charge / 1m^2 since the plate
	# has unit area.  Since we split the plate into n^2 patches, that
	# allows us to calculate dq.
	dq_top = Qtop/(n**2)
	dq_bottom = Qbottom/(n**2)

	# Calculate the field due to the top and bottom plates
	util.do_patches(x, y, z, u, v, 0.5*d, dq_top, Ex, Ey, Ez)
	util.do_patches(x, y, z, u, v, -0.5*d, dq_bottom, Ex, Ey, Ez)

	# Add physical constants.
	Ex *= 1.0/(4.0*np.pi*epsilon_0)
	Ey *= 1.0/(4.0*np.pi*epsilon_0)
	Ez *= 1.0/(4.0*np.pi*epsilon_0)

	return Ex, Ey, Ez

# Calculate the electric field on the y=0 axis as a function of x and z.
xe, ze, xc, zc = util.cartesian_2d_grid(-2, 2.0, -2, 2.0, 0.025, 0.025)
yc = np.zeros_like(xc)

Ex, Ey, Ez = charged_unit_parallelplates(xc.flatten(), yc.flatten(), zc.flatten(), -1.0, 1.0, 1.0, n=100)
Ex = np.reshape(Ex,xc.shape)
Ey = np.reshape(Ey,xc.shape)
Ez = np.reshape(Ez,xc.shape)

plt.pcolormesh(xe, ze, np.sqrt(Ex**2 + Ey**2 + Ez**2))
h=plt.colorbar()
h.set_label('|E| [V/m]')
#plt.quiver(xc,zc,Ex,Ez, color='w', scale=1e12,pivot='middle')
plt.streamplot(xc,zc,Ex,Ez, color='w')
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.xlim((-2,2))
plt.ylim((-1,1))
plt.axis('equal')
plt.show()
