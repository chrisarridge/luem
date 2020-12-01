"""Calculate the beaming pattern due to a Hertzian dipole.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches

# Don't want to rely on anyone having SciPy installed as we only want it for
# the definition of the permittivity of free space, if it's not installed then
# we just define it here.
try:
	import scipy.constants.epsilon_0 as epsilon_0
	import scipy.constants.mu_0 as mu_0
	import scipy.constants.mu_0 as c
except ImportError:
	epsilon_0 = 8.854e-12
	mu_0 = np.pi*4e-7
	c = 1/np.sqrt(epsilon_0*mu_0)

import util

# Properties of the antenna, we get p_0 from a general consideration of current
# and d: max(I)=q_0\omega, therefore, q_0 = max(I)/\omega
d = 0.01				# 1cm
Imax = 0.1				# 100 mA
omega = 2*np.pi*1e6		# 1 MHz
q_0 = Imax/omega
p_0 = d*q_0

th = np.linspace(0,2*np.pi*1.1,256)

N = lambda r, th: ((mu_0*(p_0*omega**2)**2)/(32*np.pi*np.pi*c))*(np.sin(th)/r)**2



fig = plt.figure(figsize=(12,12))
plt.polar(th, N(1.0,th))
plt.gca().set_theta_zero_location('N', offset=0)
# ax.set_rlim((1e-20,1e-10))
# ax.set_rscale('log')
plt.gca().set_theta_direction(-1)
plt.gca().set_title(r'Poynting vector [$\mathrm{W}\ \mathrm{m}^{-2}$]')
plt.savefig('hertzian_dipole.png',dpi=150)
plt.show()
