import numpy as np
import matplotlib.pyplot as pl

u0 = 10.0

xymax = 4
rmax = np.ceil(xymax*np.sqrt(2)*1.25)

fig = pl.figure(figsize=(10/2.4,10/2.4))

# Plot potential.
r, phi = np.meshgrid(np.linspace(1, rmax, 100), np.linspace(0,2*np.pi,90))
x = r*np.cos(phi)
y = r*np.sin(phi)
V = u0*(1 + r**(-2))*np.cos(phi)
h_pot_contours=pl.contourf(x, y, V, np.arange(-20,22,2))

# Plot streamlines.
r, phi = np.meshgrid(np.linspace(1, rmax, 100), np.linspace(0,2*np.pi,90))
x = r*np.cos(phi)
y = r*np.sin(phi)
psi = u0*(r - r**(-1))*np.sin(phi)
pl.contour(x, y, psi, np.concatenate((np.arange(-54,2,4),np.arange(2,58,4))), colors='w', linestyles='solid')

# Plot vector field.
r, phi = np.meshgrid(np.linspace(1, rmax, 10), np.linspace(0,2*np.pi,25))
x = r*np.cos(phi)
y = r*np.sin(phi)
ur = u0*(1-r**(-2))*np.cos(phi)
uphi = -u0*(1+r**(-2))*np.sin(phi)
ux = ur*np.cos(phi) - uphi*np.sin(phi)
uy = ur*np.sin(phi) + uphi*np.cos(phi)
pl.quiver(x,y,ux,uy)

# Finish plot.
pl.xlim([-xymax,xymax])
pl.ylim([-xymax,xymax])
pl.gca().set_aspect('equal', 'datalim')

h=pl.colorbar(h_pot_contours)
h.set_label('Potential [m m/s]')

pl.savefig('potflow.png')
pl.show()
