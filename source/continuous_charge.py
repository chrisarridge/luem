"""Plot illustration of continuous charge distributions.
"""

import numpy as np
import matplotlib.pyplot as plt
import util

# Number of point charges to generate.
n = 1000

# Total charge.
Q = 1.0

# Generate a set of randomly positioned point charges.  The positions are
# generated from a multivariate Normal distribution with a mean of x=0 and y=0,
# and a covariance of [[4,1],[1,3]].
tmp = np.random.multivariate_normal([0,0], [[4,1],[1,3]], size=n)
x = tmp[:,0]
y = tmp[:,1]

# Collect the random point charges into "continuous" rectangular patches.
dx = 0.5
dy = 0.5
xe = np.arange(-8, 8, dx)
ye = np.arange(-8, 8, dy)
xc = xe[:-1] + 0.5*dx
yc = ye[:-1] + 0.5*dy
h, _tmp1, _tmp2 = np.histogram2d(x, y, bins=(xe,ye), density=True)
h *= (Q/n)

# Plot.
fig = plt.Figure(figsize=(10/2.54,10/2.54))
plt.pcolormesh(xe, ye, h)
plt.scatter(x, y, 0.75, 'w')
plt.gca().set_aspect('equal')
# pl.gca().set_xticklabels([])
# pl.gca().set_yticklabels([])
plt.tight_layout()
# pl.savefig('continuous-charge.png',dpi=144)
plt.show()
