"""Calculate the electric field due to a finite charged rod.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import zool

# Layout a plot using the Zool library (github.com/chrisarridge/zool).
fig = zool.Figure(figwidth=zool.Fixed(20.0), figheight=zool.FromChildren(),
			layout='vertical', padding=0.5,
			margin_left=2.0, margin_right=2.0, margin_bottom=2.0, margin_top=1.0)
fig.add('base', width=zool.FromParent(), height=zool.Fixed(5.0), label='A', layout='horizontal', padding=0.5)
fig.add('base', width=zool.FromParent(), height=zool.Fixed(5.0), label='B', layout='horizontal', padding=0.5)
fig.add('A', width=zool.Fixed(14.0), height=zool.FromParent(), label='A1')
fig.add('A', width=zool.Fixed(0.5), height=zool.FromParent(), label='A2')
fig.add('B', width=zool.Fixed(14.0), height=zool.FromParent(), label='B1')
fig.add('B', width=zool.Fixed(0.5), height=zool.FromParent(), label='B2')
fig.layout()

# Lambda functions for the electric field components.
Er_fun = lambda r, z, qlambda, b: (qlambda/(4*np.pi*scipy.constants.epsilon_0)) * ((z+b)/np.sqrt(r**2 + (z+b)**2) - (z-b)/np.sqrt(r**2 + (z-b)**2))
Ez_fun = lambda r, z, qlambda, b: (qlambda/(4*np.pi*scipy.constants.epsilon_0)) * (1.0/np.sqrt(r**2 + (z-b)**2) - 1.0/np.sqrt(r**2 + (z+b)**2))

# Length and charge on the rod, then the line charge density.
b = 1.0
q = 1e8*scipy.constants.e
qlambda = q/b

# Calculate a grid.
z, r = np.meshgrid(np.linspace(-10*b,10*b,200), np.linspace(0.0,5.0,100))

# Calculate the electric field on the grid.
Er = Er_fun(r, z, qlambda, b)
Ez = Ez_fun(r, z, qlambda, b)
E = np.sqrt(Er**2 + Ez**2)

# Calculate unit field vectors.
er = Er/E
ez = Ez/E

# Make the plot.
fig.new()
ax = fig.axes('A1')
c = ax.contourf(z, r, np.log10(E), np.linspace(-3,1,100))
ax.quiver(z[::10,::10], r[::10,::10], ez[::10,::10], er[::10,::10], color='w')
ax.streamplot(z, r, Ez, Er, color=(0.75,0.75,0.75), linewidth=0.5)
h=plt.colorbar(c, fig.axes('A2'))
h.set_label(r'log$_{10}$(Electric field strength [V/m])')
h.set_ticks([-3,-2,-1,0,1])
ax.set_xlim([-10*b,10*b])
ax.set_ylim([0,5])
ax.set_ylabel('r [m]')
ax.set_xticklabels([])

ax = fig.axes('B1')
z = np.linspace(-10*b,10*b,800)
ax.semilogy(z, Er_fun(0, z, qlambda, b), label=r'$E_r$')
ax.semilogy(z, np.abs(Ez_fun(0, z, qlambda, b)), label=r'$|E_z|$')

z = np.linspace(-10*b,10*b,800)
ax.semilogy(z,0.2/z**2,'--',label=r'$z^{-2}$')
ax.legend()
ax.set_xlabel('z [m]')
ax.set_ylabel('Electric field [V/m]')
ax.set_xlim([-10*b,10*b])

#pl.savefig('electric-field.rod.png',dpi=150)
