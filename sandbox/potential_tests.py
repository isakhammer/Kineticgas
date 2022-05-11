from pykingas import KineticGas
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import Boltzmann
plt.style.use('default')

kin = KineticGas('AR,C1', potential_mode='mie')
sigma = kin.sigma_ij[0, 0]
r_list = np.linspace(sigma, 2 * sigma)
u1_list = np.empty_like(r_list)
u2_list = np.empty_like(r_list)
u12_list = np.empty_like(r_list)

F1_list = np.empty_like(r_list)
F2_list = np.empty_like(r_list)
F12_list = np.empty_like(r_list)

for i, r in enumerate(r_list):
    u1_list[i] = kin.cpp_kingas.potential(1, r, 0)
    u2_list[i] = kin.cpp_kingas.potential(2, r, 0)
    u12_list[i] = kin.cpp_kingas.potential(12, r, 0)

    F1_list[i] = kin.cpp_kingas.potential_derivative_r(1, r, 0)
    F2_list[i] = kin.cpp_kingas.potential_derivative_r(2, r, 0)
    F12_list[i] = kin.cpp_kingas.potential_derivative_r(12, r, 0)

fig, ax = plt.subplots(1, 1)
twn = ax.twinx()

plt.sca(ax)
plt.plot(r_list / sigma, u1_list / Boltzmann, label=1, color='r')
#plt.plot(r_list / sigma, u2_list / Boltzmann, label=2, color='g')
#plt.plot(r_list / sigma, u12_list / Boltzmann, label=12, color='b')
plt.ylabel('Potential')

plt.sca(twn)
plt.plot(r_list / sigma, F1_list * sigma / Boltzmann, color='r')
#plt.plot(r_list / sigma, F2_list * sigma / Boltzmann, color='g')
#plt.plot(r_list / sigma, F12_list * sigma / Boltzmann, color='b')
plt.ylabel('Force')

plt.sca(ax)
plt.legend()
plt.show()

