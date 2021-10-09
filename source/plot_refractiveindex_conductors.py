"""Plot the (complex) refractive index for conductors
"""
import numpy as np
import matplotlib.pyplot as plt

# Frequency grid.
f = np.logspace(3,12,300)

# Lambda functions for the refractive index.
n1squared = lambda f, sigma, er=1.0: 0.5*er*(1 + np.sqrt(1 + 1/((2*np.pi*f)**2 * 8.854e-12*er/sigma)) )
n2squared = lambda f, sigma, er=1.0: 0.5*er*(-1 + np.sqrt(1 + 1/((2*np.pi*f)**2 * 8.854e-12*er/sigma)) )
n1 = lambda f, sigma, er=1.0: np.sqrt(n1squared(f,sigma,er))
n2 = lambda f, sigma, er=1.0: np.sqrt(n2squared(f,sigma,er))
alpha = lambda f, sigma, er=1.0: 2.998e8/(2*np.pi*f*n2(f, sigma, er))

# Figure 20cm x 10cm (the dimensions in Matplotlib are in inches, 1"=2.54 cm).
fig = plt.figure(figsize=(20/2.54,12/2.54))

plt.subplot(2,1,1)
plt.loglog(f, n1(f, 1e-4), label=r'$10^{-4}$ S/m')
plt.loglog(f, n1(f, 1e0), label=r'$1$ S/m')
plt.loglog(f, n1(f, 1e6), label=r'$10^{6}$ S/m')
plt.ylim(1e0,1e3)
plt.plot([3e0,3e0],plt.gca().get_ylim(), '--k')
plt.plot([300e6,300e6],plt.gca().get_ylim(), '--k')
plt.plot([300e9,300e9],plt.gca().get_ylim(), '--k')
plt.plot([430e12,430e12],plt.gca().get_ylim(), '--k')
plt.plot([750e12,750e12],plt.gca().get_ylim(), '--k')
plt.plot([30e15,30e15],plt.gca().get_ylim(), '--k')
plt.plot([30e19,30e19],plt.gca().get_ylim(), '--k')
# plt.plot([3e0,300e6],plt.gca().get_ylim(), '--k')
# plt.plot([300e6,300e9],plt.gca().get_ylim(), '--k')
# plt.plot([300e9,430e12],plt.gca().get_ylim(), '--k')
# plt.plot([430e12,750e12],plt.gca().get_ylim(), '--k')
# plt.plot([750e12,30e15],plt.gca().get_ylim(), '--k')
# plt.plot([30e15,30e19],plt.gca().get_ylim(), '--k')
plt.ylabel(r'$n_1$')
plt.grid()
plt.legend()

plt.subplot(2,1,2)
plt.loglog(f, n2(f, 1e-4))
plt.loglog(f, n2(f, 1e0))
plt.loglog(f, n2(f, 1e6))
plt.ylim(1e-4,1e3)
plt.plot([3e0,3e0],plt.gca().get_ylim(), '--k')
plt.plot([300e6,300e6],plt.gca().get_ylim(), '--k')
plt.plot([300e9,300e9],plt.gca().get_ylim(), '--k')
plt.plot([430e12,430e12],plt.gca().get_ylim(), '--k')
plt.plot([750e12,750e12],plt.gca().get_ylim(), '--k')
plt.plot([30e15,30e15],plt.gca().get_ylim(), '--k')
plt.plot([30e19,30e19],plt.gca().get_ylim(), '--k')
plt.ylabel(r'$n_2$')
plt.xlabel('Frequency [Hz]')
plt.grid()
#
# plt.subplot(2,2,4)
# plt.loglog(f, alpha(f, 1e-4))
# plt.loglog(f, alpha(f, 1e0))
# plt.loglog(f, alpha(f, 1e6))
# #plt.ylim(1e-2,1e3)
# plt.plot([3e0,3e0],plt.gca().get_ylim(), '--k')
# plt.plot([300e6,300e6],plt.gca().get_ylim(), '--k')
# plt.plot([300e9,300e9],plt.gca().get_ylim(), '--k')
# plt.plot([430e12,430e12],plt.gca().get_ylim(), '--k')
# plt.plot([750e12,750e12],plt.gca().get_ylim(), '--k')
# plt.plot([30e15,30e15],plt.gca().get_ylim(), '--k')
# plt.plot([30e19,30e19],plt.gca().get_ylim(), '--k')
# plt.ylabel(r'Penetration depth, $\alpha$ [m]')
# plt.xlabel('Frequency [Hz]')
# plt.grid()
plt.savefig('refractive-index.png',dpi=150)
plt.show()