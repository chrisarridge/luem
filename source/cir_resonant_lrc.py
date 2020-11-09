"""Plot impedance for a resonant RLC circuit.
"""

import numpy as np
import matplotlib.pyplot as plt

rlc_imp = lambda omega, r, l, c: np.sqrt(r**2 + (omega*l - 1.0/(omega*c))**2)
rlc_imp_from_omega0_c = lambda omega, r, omega0, c: np.sqrt(r**2 + (1/c**2)*(omega/omega0**2 - 1.0/omega)**2)

f = np.logspace(1,5,500)
omega = 2*np.pi*f

#plt.semilogx(f, 5.0/rlc_imp_from_omega0_c(omega, 4700.0, 2*np.pi*1e3, 1e-6), label=r'R=4.7 k$\Omega$')
plt.figure(figsize=(14/2.54,12/2.54))
plt.semilogx(f, 5.0/rlc_imp_from_omega0_c(omega, 470.0, 2*np.pi*1e3, 1e-6), label=r'R=470 $\Omega$')
plt.semilogx(f, 5.0/rlc_imp_from_omega0_c(omega, 47.0, 2*np.pi*1e3, 1e-6), label=r'R=47 $\Omega$')
plt.semilogx(f, 5.0/rlc_imp_from_omega0_c(omega, 12, 2*np.pi*1e3, 1e-6), label=r'R=12 $\Omega$')
plt.xlabel(r'Driving frequency, $f$ [Hz]')
plt.ylabel(r'Current, $I$ [A]')
plt.legend()
plt.savefig('rlc-resonance.png',dpi=150)
plt.show(block=True)
