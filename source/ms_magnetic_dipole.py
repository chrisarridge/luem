import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors
import matplotlib.patches

r, th = np.meshgrid(np.linspace(1, 15, 100), np.linspace(0, 2*np.pi, 60))

V = np.cos(th)*(r**(-2))

fig = pl.figure(figsize=(16,9))
circle = matplotlib.patches.Circle((0,0), 1, edgecolor='k', facecolor='w')

ax = fig.add_subplot(141, projection='polar')
h=ax.contourf(th, r, V, 16, vmin=-1, vmax=1, cmap='RdBu_r')
ax.set_theta_zero_location('N', offset=0)
ax.set_rorigin(-1)
ax.set_theta_direction(-1)
ax.set_ylim([1,5])
ax.set_yticks([2,3,4])
h=pl.colorbar(h, orientation='horizontal')
h.set_ticks([-1,-0.5,0,0.5,1])
h.set_label('Magnetic Potential [T m]')
ax.set_title('Potential')


ax = fig.add_subplot(142, projection='polar')
h=ax.contourf(th, r, V, 16, vmin=-1, vmax=1, cmap='RdBu_r')
ax.set_theta_zero_location('N', offset=0)
ax.set_rorigin(-1)
ax.set_theta_direction(-1)
ax.set_ylim([1,5])
ax.set_yticks([2,3,4])

Br = 2*np.cos(th)*(r**(-3))
Bth = np.sin(th)*(r**(-3))
B = np.sqrt(Br**2 + Bth**2)
ax.quiver(th, r, Br*np.sin(th) + Bth*np.cos(th), Br*np.cos(th) - Bth*np.sin(th))

h=pl.colorbar(h, orientation='horizontal')
h.set_ticks([-1,-0.5,0,0.5,1])
h.set_label('Magnetic Potential [T m]')
ax.set_title('Potential with\nField Vectors')


ax = fig.add_subplot(143, projection='polar')
h=ax.contourf(th, r, V, 16, vmin=-1, vmax=1, cmap='RdBu_r')
ax.set_theta_zero_location('N', offset=0)
ax.set_rorigin(-1)
ax.set_theta_direction(-1)
ax.set_ylim([1,5])
ax.set_yticks([2,3,4])

Br = 2*np.cos(th)*(r**(-3))
Bth = np.sin(th)*(r**(-3))
B = np.sqrt(Br**2 + Bth**2)
bx = (Br*np.sin(th) + Bth*np.cos(th))/B
bz = (Br*np.cos(th) - Bth*np.sin(th))/B
ax.quiver(th[::2,::2], r[::2,::2], bx[::2,::2], bz[::2,::2], scale=25)

h=pl.colorbar(h, orientation='horizontal')
h.set_ticks([-1,-0.5,0,0.5,1])
h.set_label('Magnetic Potential [T m]')
ax.set_title('Potential with\nField Unit Vectors')



ax = fig.add_subplot(144, projection='polar')
h=ax.contourf(th, r, V, 16, vmin=-1, vmax=1, cmap='RdBu_r')
ax.set_theta_zero_location('N', offset=0)
ax.set_rorigin(-1)
ax.set_theta_direction(-1)
ax.set_ylim([1,5])
ax.set_yticks([2,3,4])

Br = 2*np.cos(th)*(r**(-3))
Bth = np.sin(th)*(r**(-3))
B = np.sqrt(Br**2 + Bth**2)
bx = (Br*np.sin(th) + Bth*np.cos(th))/B
bz = (Br*np.cos(th) - Bth*np.sin(th))/B
ax.quiver(th[::2,::2], r[::2,::2], bx[::2,::2], bz[::2,::2], scale=25)

for th_max in np.radians([10,20,30,40,50,60]):
	# th_max = np.arcsin(np.sqrt(1/r_eq))
	r_eq = 1/np.sin(th_max)**2
	th = np.arange(th_max,np.pi-th_max+np.pi/180,np.pi/180)
	pl.plot(th, r_eq*np.sin(th)**2, '-k', linewidth=2)
	pl.plot(th+np.pi, r_eq*np.sin(th)**2, '-k', linewidth=2)

h=pl.colorbar(h, orientation='horizontal')
h.set_ticks([-1,-0.5,0,0.5,1])
h.set_label('Magnetic Potential [T m]')
ax.set_title('Potential with\nField Lines')

pl.savefig('dipole.png')
pl.show(block=True)
