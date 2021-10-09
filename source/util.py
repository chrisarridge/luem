"""This module contains some utility functions that are used in other places."""

import numpy as np

# Don't want to rely on anyone having SciPy installed as we only want it for
# the definition of the permittivity of free space, if it's not installed then
# we just define it here.
try:
	import scipy.constants.epsilon_0 as epsilon_0
except ImportError:
	epsilon_0 = 8.854e-12


# Lambda function that returns the potential of a point charge.
Vmono = lambda q, x, y: q/(4*np.pi*epsilon_0*((x**2 + y**2)**(3/2)))


def cartesian_2d_grid(xmin, xmax, ymin, ymax, dx, dy):
	"""Create a 2D grid

	Each point (really a small rectangle) in a 2D grid is characterised by its
	edges in x and y (the edges of the rectangle) and centre coordintes x and y.

	Args:

		:param xmin: Minimum x value.
		:param xmax: Maximum x value.
		:param ymin: Minimum y value.
		:param ymax: Maximum y value.
		:param dx: size of each small rectangle in the x direction.
		:param dy: size of each small rectangle in the x direction.
		:return: xedges, yedges, xcentres, ycentres
	"""
	xedges = np.arange(xmin, xmax+dx, dx)
	yedges = np.arange(ymin, ymax+dy, dy)
	xe, ye = np.meshgrid(xedges, yedges)
	xc, yc = np.meshgrid((xedges[:-1]+xedges[1:])/2.0,
							(yedges[:-1]+yedges[1:])/2.0)
	return xe, ye, xc, yc



def cartesian_3d_grid(xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz):
	"""Create a 3D grid

	Each point (really a small cube) in a 3D grid is characterised by its
	edges in x, y and z (edges of the cube) and centre coordintes x, y and z.

	Args:

		:param xmin: Minimum x value.
		:param xmax: Maximum x value.
		:param ymin: Minimum y value.
		:param ymax: Maximum y value.
		:param zmin: Minimum z value.
		:param zmax: Maximum z value.
		:param dx: size of each small cube in the x direction.
		:param dy: size of each small cube in the y direction.
		:param dz: size of each small cube in the z direction.
		:return: xedges, yedges, zedges, xcentres, ycentres, zcentres
	"""
	xedges = np.arange(xmin-dz/2.0, xmax+dx/2.0, dx)
	yedges = np.arange(ymin-dz/2.0, ymax+dy/2.0, dy)
	zedges = np.arange(zmin-dz/2.0, zmax+dz/2.0, dz)
	xe, ye, ze = np.meshgrid(xedges, yedges, zedges)
	xc, yc, zc = np.meshgrid((xedges[:-1]+xedges[1:])/2.0,
				(yedges[:-1]+yedges[1:])/2.0, (zedges[:-1]+zedges[1:])/2.0,
				indexing='ij')
	return xe, ye, ze, xc, yc, zc


def do_patches(x, y, z, xp, yp, zp, dq, Ex, Ey, Ez):
	"""Adds the electric field contributions from a set of patches.

	There should be a set of N patches with centre positions (xp,yp,zp) and
	charges, then the routine calculates the relative distance for the points
	from each patch, the vector and then the field contribution.  The vector of
	charges on each patch and the positions must all be the same size, e.g.,
	if dq is 1x100, then xp, yp and zp should be 1x100.  Alternatively,
	dq can be a single value if it is constant for each patch, but then the
	patch coordinates must all be the same size. Each patch basically
	contributes a point-charge field.

	Args:

		:param x: X coordinates for the points to calculate the field at [m].
		:param y: Y coordinates for the points to calculate the field at [m].
		:param z: Z coordinates for the points to calculate the field at [m].
		:param xp: X coordinates for the centres of each patch [m].
		:param yp: Y coordinates for the centres of each patch [m].
		:param zp: Z coordinates for the centres of each patch [m].
		:param dq: Charge on each patch [C].
		:param Ex: X component of the electric field [V/m].
		:param Ey: Y component of the electric field [V/m].
		:param Ez: Z component of the electric field [V/m].
	"""
	# Loop over all the points we have been asked to calculate the electric
	# field at.
	for i in range(len(x)):
		# Cube distance from this point to each patch.
		r = np.sqrt((x[i]-xp)**2 + (y[i]-yp)**2 + (z[i]-zp)**2)

		# Electric field contribution: dq*(r-r')/|r-r'|^3 from all the patches.
		Ex[i] += np.sum(dq*(x[i]-xp)/(r**3))
		Ey[i] += np.sum(dq*(y[i]-yp)/(r**3))
		Ez[i] += np.sum(dq*(z[i]-zp)/(r**3))
