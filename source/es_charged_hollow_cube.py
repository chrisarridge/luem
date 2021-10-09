"""This code calculates the electric field due to a charged hollow cube.
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



def charged_hollow_unit_cube(x, y, z, Q, n=10):
	"""Calculates the electric field of a hollow unit cube at supplied coordinates.

	Each face of the cube is split into n^2 small rectangles (patches) and the
	total charge on the cube is divided up into these small rectangles.  Then
	the Coulomb law is used to sum up the field contributions from all the
	patches.

	Args:

		:param x: X coordinates for the points to calculate the field at [m].
		:param y: X coordinates for the points to calculate the field at [m].
		:param z: Z coordinates for the points to calculate the field at [m].
		:param Q: Total charge on the cube [C].
		:param n: Optional, number of cells to split in each dimension (default=10).
		:return: Ex [V/m], Ey [V/m], Ez [V/m]
	"""

	# Setup a grid of points over the face of a unit cube; each axis is split
	# into n chunks between -0.5 and 0.5, producing n^2 patches over the face.
	u, v = np.meshgrid(np.linspace(-0.5,0.5,n), np.linspace(-0.5,0.5,n))
	dp = u[0,1]-u[0,0]
	dS = dp*dp		# area of the patch - these are all uniform.
	u = u.flatten()
	v = v.flatten()

	# Clear the electric field to start the summation.
	Ex = np.zeros(len(x))
	Ey = np.zeros(len(x))
	Ez = np.zeros(len(x))

	# The surface charge density is just 1/6th of the total charge since each
	# face has unit area.  Since we split each face into n^2 patches, that
	# allows us to calculate dq.
	dq = (Q/6)/(n**2)

	# top face (+z) and bottom face (-z)
	util.do_patches(x, y, z, u, v, 0.5, dq, Ex, Ey, Ez)
	util.do_patches(x, y, z, u, v, -0.5, dq, Ex, Ey, Ez)

	# front (+x) and back face (-x)
	util.do_patches(x, y, z, 0.5, u, v, dq, Ex, Ey, Ez)
	util.do_patches(x, y, z, -0.5, u, v, dq, Ex, Ey, Ez)

	# left (-y) and right face (+y)
	util.do_patches(x, y, z, u, 0.5, v, dq, Ex, Ey, Ez)
	util.do_patches(x, y, z, u, -0.5, v, dq, Ex, Ey, Ez)

	# Add physical constants.
	Ex *= 1.0/(4.0*np.pi*epsilon_0)
	Ey *= 1.0/(4.0*np.pi*epsilon_0)
	Ez *= 1.0/(4.0*np.pi*epsilon_0)

	return Ex, Ey, Ez

# Calculate the electric field on the z=0 axis as a function of x and y.
xe, ye, xc, yc = util.cartesian_2d_grid(-0.6, 0.6, -0.6, 0.6, 0.005, 0.005)
zc = np.zeros_like(xc)

Ex, Ey, Ez = charged_hollow_unit_cube(xc.flatten(), yc.flatten(), zc.flatten(), 1, n=20)
Ex = np.reshape(Ex,xc.shape)
Ey = np.reshape(Ey,xc.shape)
Ez = np.reshape(Ez,xc.shape)

#plt.pcolormesh(xe, ye, np.sqrt(Ex**2 + Ey**2))
#plt.colorbar()
#plt.pcolormesh(xe, ye, np.sqrt(Ex**2 + Ey**2))
#plt.colorbar()
#plt.clim([0,1e9])
#plt.quiver(xc,yc,Ex,Ey)
plt.streamplot(xc,yc,Ex,Ey,density=2)
plt.plot([-0.5,0.5],[-0.5,-0.5],'-k')
plt.plot([-0.5,0.5],[0.5,0.5],'-k')
plt.plot([-0.5,-0.5],[-0.5,0.5],'-k')
plt.plot([0.5,0.5],[-0.5,0.5],'-k')
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.gca().set_aspect('equal')
plt.axis('equal')
plt.show()
