"""Plot the potential due to a charged infinite sheet.
"""

import numpy as np
import matplotlib.pyplot as plt

z = np.linspace(-10,10,200)

fig = plt.figure(figsize=(20/2.54,10/2.54))

# Plot the potential.
plt.subplot(2,1,1)
plt.plot(z, -0.5*np.abs(z))
plt.gca().set_xticklabels([])
plt.grid()
plt.ylabel(r'Potential, V(z) [$\frac{\sigma}{\epsilon_0}$ V]')

# Plot the electric field.
plt.subplot(2,1,2)
plt.plot(z, 0.5*np.sign(z))
plt.grid()
plt.ylabel(r'Electric field,'+'\n'+r'$\mathrm{E}_\mathrm{z}$(z) [$\frac{\sigma}{\epsilon_0}$ V/m]')
plt.xlabel('z [m]')
plt.show()
