"""Plot the velocity field around a Rankine vortex and the curl.
"""
import numpy as np
import matplotlib.pyplot as plt
import util

def rankine_uph(r, k, a):
	"""Return the azimuthal part of the velocity field for a Rankine vortex.

	Args:

		:param r: Radial coordinate.
		:param k: Constant.
		:param a: Radius of the vortex.
		:returns: phi component of the velocity field.
	"""
	uph = k*r/a
	uph[r>a] = k*a/r[r>a]
	return uph

def rankine_curl(r, k, a):
	"""Return the z component of the curl of a Rankine vortex velocity field.

	Args:

		:param r: Radial coordinate.
		:param k: Constant.
		:param a: Radius of the vortex.
		:returns: z component of the curl of the velocity field.
	"""
	curl_in_z = np.zeros_like(r)
	curl_in_z[r<=a] = 2*k*a
	return curl_in_z



# Create a 2D grid in Cartesian and calculate the position in cylindrical polar.
xe, ye, xc, yc = util.cartesian_2d_grid(-4.5, 4.5, -3.5, 3.5, 0.02, 0.02)
r = np.sqrt(xc**2 + yc**2)
ph = np.arctan2(yc, xc)

# Plot the curl.
curl = rankine_curl(r, 1.0, 1.0)
plt.pcolormesh(xe, ye, curl, zorder=-1, cmap='RdBu')
h=plt.colorbar()
h.set_label(r'$(\nabla\times\mathbf{u})_z$')

# Calculate the velocity field and convert the vectors back to Cartesian.
uph = rankine_uph(r, 1.0, 1.0)
ux = -uph*np.sin(ph)
uy = uph*np.cos(ph)

# Plot the stream (field) lines and vectors.
plt.streamplot(xc, yc, ux, uy, color=(0.9,0.9,0.9), zorder=0)
plt.quiver(xc[::20,::20], yc[::20,::20], ux[::20,::20], uy[::20,::20], pivot='middle', zorder=1)

plt.xlim([-4,4])
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal','datalim')
#plt.savefig('rankine.png',dpi=150)
plt.show()
