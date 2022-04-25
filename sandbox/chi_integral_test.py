from pykingas import KineticGas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from scipy.constants import Boltzmann
plt.style.use('default')

T = 300
kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]

b_list = np.linspace(0, 15, 30) * sigma
g_list = np.linspace(1e-5, 5, 30)

chi_list = np.empty((len(b_list), len(g_list)))

f = lambda x, y: x + y**2

for bi, b in enumerate(b_list):
    for gi, g in enumerate(g_list):
        print(b/sigma, g)
        chi_list[gi, bi] = kin.cpp_kingas.chi(1, T, g, b)

b_list, g_list = np.meshgrid(b_list, g_list)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_wireframe(b_list/sigma, g_list, chi_list / np.pi)
ax.set_xlabel(r'$b$ [$\sigma$]')
ax.set_ylabel(r'$g$ [-]')
ax.set_zlabel(r'$\chi$ [$\pi$]')
plt.show()
