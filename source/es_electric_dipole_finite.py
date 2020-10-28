"""Calculate the electric field due to a dipole but where each charge
was a finite size.
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


def dipole_rho(x, y, x1, y1, x2, y2, qrho1, qrho2, a1, a2):
	"""Calculate the charge density due to a pair of uniformly charged spheres

	Args:

		:param x: X coordinate(s) to calculate the charge density at.
		:param y: Y coordinate(s) to calculate the charge density at.
		:param x1: X coordinate of the centre of the first sphere.
		:param y1: Y coordinate of the centre of the first sphere.
		:param x2: X coordinate of the centre of the second sphere.
		:param y2: Y coordinate of the centre of the second sphere.
		:param qrho1: (Uniform) charge density of the first sphere.
		:param qrho2: (Uniform) charge density of the second sphere.
		:param a1: Radius of the first sphere.
		:param a2: Radius of the second sphere.
		:return: Charge density at each (x,y) point.
	"""
	r1 = np.sqrt((x-x1)**2 + (y-y1)**2)
	r2 = np.sqrt((x-x2)**2 + (y-y2)**2)
	rho = np.zeros_like(r1)
	rho[r1<=a1] = qrho1
	rho[r2<=a2] = qrho2
	return rho


def dipole_field(x, y, x1, y1, x2, y2, qrho1, qrho2, a1, a2):
	"""Calculate the electric field due to a pair of uniformly charged spheres

	Args:

		:param x: X coordinate(s) to calculate the charge density at.
		:param y: Y coordinate(s) to calculate the charge density at.
		:param x1: X coordinate of the centre of the first sphere.
		:param y1: Y coordinate of the centre of the first sphere.
		:param x2: X coordinate of the centre of the second sphere.
		:param y2: Y coordinate of the centre of the second sphere.
		:param qrho1: (Uniform) charge density of the first sphere.
		:param qrho2: (Uniform) charge density of the second sphere.
		:param a1: Radius of the first sphere.
		:param a2: Radius of the second sphere.
		:return: Electric field in the x and y directions at each (x,y) point.
	"""
	# Convert the cartesian coordinates into polar coordinates centred on each
	# charged sphere.
	r1 = np.sqrt((x-x1)**2 + (y-y1)**2)
	ph1 = np.arctan2(y-y1, x-x1)
	r2 = np.sqrt((x-x2)**2 + (y-y2)**2)
	ph2 = np.arctan2(y-y2, x-x2)

	# Calculate the total charge on each sphere from the charge density, the
	# radius and hence the volume.
	q1 = qrho1*(4/3)*np.pi*a1**3
	q2 = qrho2*(4/3)*np.pi*a2**3

	# Calculate the electric field due to each sphere (note this is radially
	# out from each sphere.
	Er1 = (q1/(4*np.pi*epsilon_0))*r1/(a1**3)
	Er1[r1>a1] = (q1/(4*np.pi*epsilon_0))/r1[r1>a1]**2
	Er2 = (q2/(4*np.pi*epsilon_0))*r2/(a2**3)
	Er2[r2>a2] = (q2/(4*np.pi*epsilon_0))/r2[r2>a2]**2

	# Convert to Cartesian coordinates and return the vectors.
	return Er1*np.cos(ph1) + Er2*np.cos(ph2), Er1*np.sin(ph1) + Er2*np.sin(ph2)



# Side dimension [m] of the dipole and magnitude [C] of the charges.
d = 5.0
qmag = 1.0
a = 1.0



# Setup a plot that is 16 cm x 12 cm (sizes need to be give in inches).
fig = plt.figure(figsize=(16/2.54,12/2.54))

# Calculate the charge density on the grid and plot this.
xe, ye, xc, yc = util.cartesian_2d_grid(-10, 10, -10, 10, 0.025, 0.025)
rho = dipole_rho(xc, yc, -d, 0, d, 0, qmag, -qmag, a, a)
plt.pcolormesh(xe, ye, rho, cmap='RdBu', zorder=-10)
h=plt.colorbar()
h.set_label(r'Charge density (also proportional to div E)[$\mathrm{C}\ \mathrm{m}^{-3}$]')

# Calculate the electric field and unit field.
Ex, Ey = dipole_field(xc, yc, -d, 0, d, 0, qmag, -qmag, a, a)
ex, ey = Ex/np.sqrt(Ex**2+Ey**2), Ey/np.sqrt(Ex**2+Ey**2)

# Plot field vectors; we slice the arrays to reduce the number of vectors
plt.quiver(xc[::19,::19], yc[::19,::19], Ex[::19,::19], Ey[::19,::19], pivot='middle', zorder=0)

# Plot field lines.
plt.streamplot(xc, yc, Ex, Ey, color=(0.5,0.5,0.5), zorder=-5)

# Finish annotations.
plt.axis('equal')
# plt.tight_layout()
plt.xlim([-10,10])
# plt.xlim([-10,10])
plt.xlabel('x [m]')
plt.ylabel('x [m]')

# plt.savefig('electric_dipole.png',dpi=144)
plt.show(block=True)


plt.show()
