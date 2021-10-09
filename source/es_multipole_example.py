"""This code sets up a dipole, monopole and quadrupole and takes moments
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

# Define the monopole in terms of the charge (q) and the "size" of the charge clouds.
q_monopole = 10
a_monopole = 1.0

# Define the dipole in terms of the charge (q), the distance between charges (d)
# and the size of the charge clouds, also calculate the analytical value of the
# dipole moment.
q_dipole = 1.0
d_dipole = 1.0
p_dipole = q_dipole*d_dipole
a_dipole = 0.4

# Define the quadrupole in terms of the charge (q), the distance between charges (d)
# and the size of the charge clouds.
q_quad = 0.5
a_quad = 0.2
d_quad = 3.0


# Each "charge" is modelled using a 3D Gaussian:
# f(x,y,z) = q*(a*sqrt(2*pi)^-3) exp[-(0.5/a^2)*((x-x0)^2-(y-y0)^2-(z-z0)^2)]
#
# where (x0,y0,z0) is the position of the charge, q is the total charge, and a
# is a scale defining the "width" of the charge in space - larger values mean
# the charge is more spread out.
f = lambda x,y,z, x0, y0, z0, q, a: q * (np.sqrt(2*np.pi)*a)**(-3.0) * np.exp(-(0.5/a**2)*((x-x0)**2+(y-y0)**2+(z-z0)**2))

# Setup grid of points to calculate charge density.
dx, dy, dz = 0.05, 0.05, 0.05
xe, ye, ze, xc, yc, zc = util.cartesian_3d_grid(-5,5, -5,5, -5,5, dx, dy, dz)

# Calculate the charge density due to the monopole, the dipole and the
# quadrupole and sum them together.
rho_monopole = f(xc, yc, zc, 0.0, 0.0, 0.0, q_monopole, a_monopole)
rho_dipole = (f(xc, yc, zc, 0.0, -d_dipole/2.0, 0.0, -q_dipole, a_dipole) +
				f(xc, yc, zc, 0.0, d_dipole/2.0, 0.0, q_dipole, a_dipole))
rho_quad = (f(xc, yc, zc, -d_quad/2.0, -d_quad/2.0, 0.0, -q_quad, a_quad) +
			f(xc, yc, zc, -d_quad/2.0, d_quad/2.0, 0.0, q_quad, a_quad) +
			f(xc, yc, zc, d_quad/2.0, -d_quad/2.0, 0.0, q_quad, a_quad) +
			f(xc, yc, zc, d_quad/2.0, d_quad/2.0, 0.0, -q_quad, a_quad))
rho = rho_monopole + rho_dipole + rho_quad


# ==============================================================================
#
# Plot the charge densities.
#
# ==============================================================================
fig = plt.figure(figsize=(16.0,4.0))

vminmax = 1

plt.subplot(1,4,1)
plt.title('Monopole')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], rho_monopole[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')
plt.ylabel('y [m]')

plt.subplot(1,4,2)
plt.title('Dipole')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], rho_dipole[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

plt.subplot(1,4,3)
plt.title('Quadrupole')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], rho_quad[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

plt.subplot(1,4,4)
plt.title('Total')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], rho[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

plt.tight_layout()
fig.subplots_adjust(right=0.8)#

pos = plt.gca().get_position()
cbar_ax = fig.add_axes([0.85, pos.y0, 0.02, pos.height])
h=plt.colorbar(cax=cbar_ax)
h.set_label(r'Charge density, $\rho\ [\rm{C}\ \rm{m}^{-3}]$')

plt.savefig('mp-example-charge-density.png',dpi=300)
plt.show()





monopole_calc = np.sum(rho)*dx*dy*dz
print('Monopole moment: {} actual={}'.format(monopole_calc,q_monopole))


# ==============================================================================
#
# Plot the dipole moment integrand.
#
# ==============================================================================
fig = plt.figure(figsize=(16.0,4.0))

vminmax = 1

plt.subplot(1,3,1)
plt.title(r'$\rho({\bf r}\prime)r_x\prime$')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], xc[:,:,100].T*rho[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')
plt.ylabel('y [m]')

plt.subplot(1,3,2)
plt.title(r'$\rho({\bf r}\prime)r_y\prime$')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], yc[:,:,100].T*rho[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

plt.subplot(1,3,3)
plt.title(r'$\rho({\bf r}\prime)r_z\prime$')
plt.pcolormesh(xe[:,:,100], ye[:,:,100], zc[:,:,100].T*rho[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

plt.tight_layout()
fig.subplots_adjust(right=0.8)#
pos = plt.gca().get_position()

cbar_ax = fig.add_axes([0.85, pos.y0, 0.02, pos.height])
h=plt.colorbar(cax=cbar_ax)
h.set_label(r'Dipole moment integrand $[\rm{C}\ \rm{m}^{-2}]$')

plt.savefig('mp-example-dipole-integrand.png',dpi=300)

plt.show()

dipole_calc_x = np.sum(xc*rho)*dx*dy*dz
dipole_calc_y = np.sum(yc*rho)*dx*dy*dz
dipole_calc_z = np.sum(zc*rho)*dx*dy*dz
print('Dipole moment: ({},{},{}) actual=({},{},{})'.format(dipole_calc_x,dipole_calc_y,dipole_calc_z,0,p_dipole,0))


# ==============================================================================
#
# Plot the quadrupole moment integrand.
#
# ==============================================================================
fig = plt.figure(figsize=(16.0,14.0))

vminmax = 1

r2 = xc**2 + yc**2 + zc**2

plt.subplot(3,3,1)
plt.title(r'$\rho({\bf r\prime})(3r_x\prime r_x\prime - (r\prime)^2)$')
ig = (3*xc*xc - r2*1.0)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.gca().set_aspect('equal')
plt.ylabel('y [m]')

plt.subplot(3,3,2)
plt.title(r'$\rho({\bf r\prime})(3r_x\prime r_y\prime)$')
ig = (3*xc*yc)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.gca().set_aspect('equal')

plt.subplot(3,3,3)
plt.title(r'$\rho({\bf r\prime})(3r_x\prime r_z\prime)$')
ig = (3*xc*zc)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.gca().set_aspect('equal')
ax3 = plt.gca()

plt.subplot(3,3,4)
plt.title(r'$\rho({\bf r\prime})(3r_y\prime r_x\prime)$')
ig = (3*yc*xc)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.gca().set_aspect('equal')
plt.ylabel('y [m]')

plt.subplot(3,3,5)
plt.title(r'$\rho({\bf r\prime})(3r_y\prime r_y\prime - (r\prime)^2)$')
ig = (3*yc*yc - r2)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.gca().set_aspect('equal')

plt.subplot(3,3,6)
plt.title(r'$\rho({\bf r\prime})(3r_y\prime r_z\prime)$')
ig = (3*yc*zc)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.gca().set_aspect('equal')

plt.subplot(3,3,7)
plt.title(r'$\rho({\bf r\prime})(3r_z\prime r_x\prime)$')
ig = (3*zc*xc)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')
plt.ylabel('y [m]')

plt.subplot(3,3,8)
plt.title(r'$\rho({\bf r\prime})(3r_z\prime r_y\prime)$')
ig = (3*zc*yc)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

plt.subplot(3,3,9)
plt.title(r'$\rho({\bf r\prime})(3r_z\prime r_z\prime - (r\prime)^2)$')
ig = (3*zc*zc - r2)*rho
plt.pcolormesh(xe[:,:,100], ye[:,:,100], ig[:,:,100].T, cmap='RdBu_r', vmin=-vminmax, vmax=vminmax)
plt.xlabel('x [m]')
plt.gca().set_aspect('equal')

ax9 = plt.gca()
plt.tight_layout()

fig.subplots_adjust(right=0.8)#
pos3=ax3.get_position()
pos9=ax9.get_position()

cbar_ax = fig.add_axes([0.85, pos9.y0, 0.02, pos3.y1-pos9.y0])
h=plt.colorbar(cax=cbar_ax)
h.set_label(r'Quadrupole moment integrand $[\rm{C}\ \rm{m}^{-1}]$')


plt.savefig('mp-example-quad-integrand.png',dpi=300)

plt.show()

quad_calc_xx = np.sum((3*xc*xc - r2)*rho)*dx*dy*dz
quad_calc_xy = np.sum((3*xc*yc)*rho)*dx*dy*dz
quad_calc_xz = np.sum((3*xc*zc)*rho)*dx*dy*dz

quad_calc_yx = np.sum((3*yc*xc)*rho)*dx*dy*dz
quad_calc_yy = np.sum((3*yc*yc - r2)*rho)*dx*dy*dz
quad_calc_yz = np.sum((3*yc*zc)*rho)*dx*dy*dz

quad_calc_zx = np.sum((3*zc*xc)*rho)*dx*dy*dz
quad_calc_zy = np.sum((3*zc*yc)*rho)*dx*dy*dz
quad_calc_zz = np.sum((3*zc*zc - r2)*rho)*dx*dy*dz
print('Quadrupole moment tensor')
print('xx={} xy={} xz={}'.format(quad_calc_xx,quad_calc_xy,quad_calc_xz))
print('yx={} yy={} yz={}'.format(quad_calc_yx,quad_calc_yy,quad_calc_yz))
print('zx={} zy={} zz={}'.format(quad_calc_zx,quad_calc_zy,quad_calc_zz))
