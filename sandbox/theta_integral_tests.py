from pykingas import KineticGas, logspace, erfspace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize, LogNorm
from scipy.integrate import quad
from scipy.special import erf
plt.style.use('default')

kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]

T, g, b = 300, 2, 7.5 * sigma

func = lambda r: kin.cpp_kingas.theta_integrand(1, T, r, g, b)
integral = lambda R_min, r: quad(func, R_min, r)[0]
theta = lambda R_min, N: kin.cpp_kingas.theta(1, T, R_min, g, b, N)

R = kin.cpp_kingas.get_R(1, T, g, b)
R_real = R / (1 + 1e-5)
r_max = R
while func(r_max) > 1e-6 * func(R):
    r_max += R

def py_logspace(lmin, lmax, N_gridpoints):
    grid = np.empty(N_gridpoints)
    dx = (lmax - lmin) / (N_gridpoints - 1)
    A = (lmin - lmax) / np.log(lmax / lmin)
    B = lmin - A * np.log(lmax)
    linearmap = lambda x: A * x + B
    for i in range(N_gridpoints):
        x = np.log(lmax - dx * i)
        grid[i] = linearmap(x)

    return grid

def py_gridspace(lmin, lmax, N_gridpoints, f):
    grid = np.empty(N_gridpoints)
    dx = (lmax - lmin) / (N_gridpoints - 1)
    A = (lmin - lmax) / (f(lmax) - f(lmin))
    B = lmin - A * f(lmax)
    for i in range(N_gridpoints):
        x = lmax - dx * i
        grid[i] = A * f(x) + B

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

#fig, axs = plt.subplots(1, 1)
#r_list = np.linspace(R, r_max)
#plt.plot(r_list, (1 + 1e-3) * r_max * erf(2 * (r_list - R) / (r_max - R)))
#plt.show()
#exit(0)

def plot_erfspace():
    N = 10
    linpoints = np.linspace(R, r_max, N)
    expo = 1
    prefac = 2
    f = lambda x: erf(prefac * ((x**expo) - (R ** expo)) / ((r_max**expo) - (R ** expo)))
    gridpoints = py_gridspace(R, r_max, N, f)
    plt.plot(linpoints, gridpoints, label=r'$A =$ ' + str(round(prefac, 2)) + r', $B =$' + str(round(expo, 2)))
    while (abs(func(gridpoints[0]) - abs(func(gridpoints[1]))) / func(gridpoints[0]) > 0.1):
        expo *= 0.9
        prefac *= 1.5
        f = lambda x: erf(prefac * ((x**expo) - (R ** expo)) / ((r_max**expo) - (R ** expo)))
        gridpoints = py_gridspace(R, r_max, N, f)
        plt.plot(linpoints, gridpoints, label=r'$A =$ '+ str(round(prefac, 2))+r', $B =$'+str(round(expo, 2)))

    print(prefac, expo)
    erfpoints = erfspace(R, r_max, N, prefac, expo)
    print(gridpoints - erfpoints)
    plt.plot(linpoints, erfpoints, label=r'$A =$ '+ str(round(prefac, 2))+r', $B =$'+str(round(expo, 2)) , marker='o')
    plt.legend()
    plt.show()


def plot_theta():
    N_list = [10, 50, 100, 200]
    markers = ['o', '.', 'x', '+', 'v', '^']
    fig, axs = plt.subplots(2, 1, sharex='all')
    ax1, ax2 = axs
    cmap = get_cmap('cool')
    norm = Normalize(vmin=min(N_list), vmax=max(N_list))
    r_min_plot = R # 7.5 * sigma
    r_max = 100 * sigma
    for i, N in enumerate(N_list):
        print('N =', N)
        expo = 1
        prefac = 2
        gridpoints = erfspace(r_min_plot, r_max, N, prefac, expo)
        while abs(func(gridpoints[0]) - abs(func(gridpoints[1]))) / func(gridpoints[0]) > 0.1:
            expo *= 0.9
            prefac *= 1.5
            gridpoints = erfspace(r_min_plot, r_max, N, prefac, expo)

        dtdr = [func(r) / np.pi for r in gridpoints]
        print('dtdr : ', dtdr)
        ax1.plot(gridpoints / sigma, dtdr, color=cmap(norm(N)), linestyle='', marker='.', label=N)

        t = theta(r_min_plot, N) # N is automatically adjusted inside theta, that is why theta(R, 10) gives the same result as theta(R, 100)
        #print()
        #print('A, B =', prefac, expo)
        #print('r_max =', r_max)
        #print('R =', R)
        print('theta0 =', func(r_min_plot))
        #print()
        ax2.plot(gridpoints / sigma, [t / np.pi for _ in gridpoints], color=cmap(norm(N)), linestyle='--')
        ax2.plot(gridpoints[2:] / sigma, [py_trapz(gridpoints, i) / np.pi for i in range(2, len(gridpoints))], color=cmap(norm(N)))

    ax1.legend()
    ax2.set_xscale('log')
    ax1.set_ylabel(r'd$\theta$/d$r$ [m$^{-1}$]')
    ax2.set_ylabel(r'$\theta$ [$\pi$]')
    ax1.set_yscale('log')
    plt.show()


def plot_theta_av_b(b_min, b_max):
    b_list = np.linspace(b_min, b_max, 50) * sigma
    R_factors = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
    norm = LogNorm(vmin=min(R_factors), vmax=max(R_factors))
    cmap = get_cmap('cool')
    for Rf in R_factors:
        t_list = np.empty_like(b_list)
        chi_list = np.empty_like(b_list)
        for i, bi in enumerate(b_list):
            R = kin.cpp_kingas.get_R(1, T, g, bi) * (1 + Rf) / (1 + 1e-5)
            t_list[i] = kin.cpp_kingas.theta(1, T, R, g, bi, 50)
        for i in range(len(b_list)):
            chi_list[i] = kin.cpp_kingas.chi(1, T, g, b_list[i]) # np.pi - 2 * (t_list[i] - t_list[-1] + np.pi / 2)#
        plt.plot(b_list / sigma, t_list / np.pi, color=cmap(norm(Rf)), label=r'$\theta$')
        plt.plot(b_list / sigma, chi_list / np.pi, color=cmap(norm(Rf)), label=r'$\chi$')
    plt.plot(b_list / sigma, np.ones_like(b_list) * 0.5, linestyle='--', color='black')
    plt.plot(b_list / sigma, np.zeros_like(b_list), linestyle='--', color='black')
    plt.legend()
    plt.xlabel(r'$b$ [$\sigma$]')
    plt.ylabel(r'$\theta$ [$\pi$]')
    plt.show()

def plot_theta_integrand_dblderiv():
    r_grid = np.array(erfspace(R, 10 * sigma, 50, 2, 1))
    integrand = np.empty_like(r_grid)
    for i, r in enumerate(r_grid):
        integrand[i] = func(r)

    d2tdr2 = np.empty(len(integrand) - 2)
    for i in range(1, len(integrand) - 1):
        t_1 = func(r_grid[i - 1])
        t0 = func(r_grid[i])
        t1 = func(r_grid[i + 1])
        h1 = r_grid[i + 1] - r_grid[i]
        h2 = r_grid[i] - r_grid[i - 1]
        d2tdr2[i - 1] = 2 * (t1 + (h1 / h2) * t_1 - (h1 / h2 + 1) * t0) / (np.power(h1, 2) + h1 * h2)

    plt.plot(r_grid[1 : -1] / sigma, d2tdr2)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

plot_theta_integrand_dblderiv()
#plot_theta_av_b(1, 15)
#plot_theta()