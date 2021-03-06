import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from coordinate_generator import radial_coords, zed_coords

r_OC, r_YMC, r_GC = radial_coords()
z_OC, z_YMC, z_GC = zed_coords()

'''
y_r_OC, x_r_OC = np.histogram(r_OC, bins=0, density=True)
y_r_YMC, x_r_YMC = np.histogram(r_YMC, bins=4, density=True)
y_r_GC, x_r_GC = np.histogram(r_GC, bins=40, density=True)

y_z_OC, x_z_OC = np.histogram(z_OC, bins=20, density=True)
y_z_YMC, x_z_YMC = np.histogram(z_YMC, bins=20, density=True)
y_z_GC, x_z_GC = np.histogram(z_GC, bins=20, density=True)

N_OC, N_YMC, N_GC = len(x_r_OC), len(x_r_YMC), len(x_r_GC)
'''

plt.figure()
plt.scatter(r_OC, np.abs(z_OC), s=1, label='Open Clusters')
plt.scatter(r_YMC, np.abs(z_YMC), s=1, label='Young Massive Clusters')
plt.scatter(r_GC, np.abs(z_GC), s=1, label='Globular Clusters')
plt.legend(loc='best')
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.xlabel(r'$r$ (kpc)', fontsize=16)
plt.ylabel(r'$|z|$ (kpc)', fontsize=16)
plt.title('Star Cluster Distribution', fontsize=12)
plt.tight_layout()
plt.savefig('clusters_scatter.pdf')
plt.close()


'''
plt.figure()
plt.plot( [ 0.5*(x_r_OC[i-1]+x_r_OC[i]) for i in range(1,N_OC) ], y_r_OC, label='Open Clusters')
plt.plot( [ 0.5*(x_r_YMC[i-1]+x_r_YMC[i]) for i in range(1,N_YMC) ], y_r_YMC, label='Young Massive Clusters')
plt.plot( [ 0.5*(x_r_GC[i-1]+x_r_GC[i]) for i in range(1,N_GC) ], y_r_GC, label='Globular Clusters')
plt.xlabel(r'$r$ (kpc)', fontsize=16)
plt.ylabel(r'Probability Density', fontsize=16)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('radial_coord_distribution.pdf')
plt.close()

plt.figure()
plt.plot( [ 0.5*(x_z_OC[i-1]+x_z_OC[i]) for i in range(1,N_OC) ], y_z_OC, label='Open Clusters')
plt.plot( [ 0.5*(x_z_YMC[i-1]+x_z_YMC[i]) for i in range(1,N_YMC) ], y_z_YMC, label='Young Massive Clusters')
plt.plot( [ 0.5*(x_z_GC[i-1]+x_z_GC[i]) for i in range(1,N_GC) ], y_z_GC, label='Globular Clusters')
plt.xlabel(r'$z$ (kpc)', fontsize=16)
plt.ylabel(r'Probability Density', fontsize=16)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('zed_coord_distribution.pdf')
plt.close()
'''
