from pykingas import KineticGas
from pykingas.KineticGas import logspace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from scipy.integrate import quad

kin = KineticGas('AR,C1', potential_mode='mie')
sigma = kin.sigma_ij[0, 0]

T, g, b = 300, 1, 1.2 * sigma

func = lambda r: kin.cpp_kingas.theta_integrand(1, T, r, g, b)
integral = lambda R_min, r: quad(func, R_min, r)[0]
theta = lambda R_min, r: kin.cpp_kingas.theta(1, T, R_min, g, b)

R = kin.cpp_kingas.get_R(1, T, g, b)
r_max = R
while func(r_max) > 1e-3 * func(R):
    r_max += R


def py_logspace(lmin, lmax, N_gridpoints):
    grid = np.empty(N_gridpoints)
    dx = (lmax - lmin) / N_gridpoints
    A = (lmin - lmax) / np.log(lmax / lmin)
    B = lmin - A * np.log(lmax)
    linearmap = lambda x: A * x + B
    for i in range(N_gridpoints):
        x = np.log(lmax - dx * i)
        grid[i] = linearmap(x)

    return grid

def py_trapz(r_grid, N_points):
    val = 0
    r1 = r_grid[0]
    r2 = r_grid[1]
    integrand1 = func(r1)
    integrand2 = func(r2)
    A_coeff = (integrand2 - integrand1) / (r2 - r1)
    B_coeff = integrand1 - A_coeff * r1
    val += A_coeff * (np.power(r2, 2) - np.power(r1, 2)) / 2 + B_coeff * (r2 - r1)
    for i in range(2, N_points):
        r1 = r_grid[i - 1]
        r2 = r_grid[i]
        integrand1 = integrand2
        integrand2 = func(r2)
        A_coeff = (integrand2 - integrand1) / (r2 - r1)
        B_coeff = integrand1 - A_coeff * r1
        val += A_coeff * (np.power(r2, 2) - np.power(r1, 2)) / 2 + B_coeff * (r2 - r1)
    return val


N_list = [100, 500, 1000]

cmap = get_cmap('cool')
norm = Normalize(vmin=min(N_list), vmax=max(N_list))
for N in N_list:
    r_list = np.logspace(np.log(R), np.log(r_max), 1000, base=np.e)
    integration_points = logspace(R, r_max, N)

    plt.plot(integration_points, [theta(R, r) for r in integration_points], color=cmap(norm(N)), linestyle='--')
    plt.plot(integration_points[1:], [py_trapz(integration_points, i) for i in range(1, N)], color=cmap(norm(N)))
    #ax2.set_xscale('log')

plt.show()