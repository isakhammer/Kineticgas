from pykingas import KineticGas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from scipy.constants import Boltzmann, Avogadro
import warnings

def xy_to_rt(x, y):
    r = np.sqrt(x**2 + y**2)
    t = np.arccos(x / r)
    return r, t

def rt_to_xy(r, t):
    x = r * np.cos(t)
    y = r * np.sin(t)
    return x, y

def vec_len(vec):
    return np.sqrt(np.sum(vec**2))

def normalize_vec(vec):
    return vec / vec_len(vec)

def printarr(arr):
    for line in arr:
        for x in line:
            print(x, end=' '*(20 - len(str(x))))
        print()

def total_energy(r, t, g):
    return kin.cpp_kingas.potential(1, r, t) * kin.m0 / (np.prod(kin.mole_weights)) + 0.5 * vec_len(g)**2

def potential_energy(r, t):
    return kin.cpp_kingas.potential(1, r, t) * kin.m0 / (np.prod(kin.mole_weights))

def get_path(T, g, b, y0=5):
    g = g * np.sqrt(2 * Boltzmann * T * kin.m0 / np.prod(kin.mole_weights))
    print('g_real =', round(g, 2), 'm/s')
    sigma = kin.sigma_ij[0, 0]
    y0 = y0 * sigma
    b = b * sigma
    g = np.array([0, - g])  # Rett nedover
    x = b
    y = y0
    r0, t = xy_to_rt(x, y)
    r = r0

    x_list = [x]
    y_list = [y]
    g_list = [vec_len(g)]
    E_list = [total_energy(r, t, g)]

    F = kin.cpp_kingas.potential_derivative_r(1, r, t) * kin.m0 / (np.prod(kin.mole_weights))
    F_vec = - F * normalize_vec(np.array([x, y]))
    a = F_vec

    dt = - 0.1 * (sigma / g[1]) # 10 tidssteg for å bevege seg 1 sigma
    i = 0
    while r <= r0:
        pos = np.array([x_list[i], y_list[i]])  # Posisjon
        r, t = xy_to_rt(pos[0], pos[1])
        if (np.dot(g, normalize_vec(pos)) < 0 # Partikkelen er på vei mot potensialet
                and (E_list[0] - potential_energy(r, t) < 0.05 * abs(E_list[0])) # Potensiell energi er veldig stor
                and vec_len(g) * dt < 5e-2 * sigma): # Tidssteg er veldig lite
            print('HS at', pos)
            g = g - 2 * normalize_vec(pos) * np.dot(g, normalize_vec(pos)) # Behandle potensialet som en hard kule (speil g-vektor om pos-vektor)


        pos += g * dt  # Ny posisjon
        r, t = xy_to_rt(pos[0], pos[1])

        if potential_energy(r, t) > E_list[0]: #Sørger for energibevaring
            dt *= 0.5 # Reduser tidssteg og beregn forflytning på nytt
        else:
            g = g + a * dt  # Ny hastighet
            g = normalize_vec(g) * np.sqrt(2 * (E_list[0] - potential_energy(r, t))) # Korrigerer for energibevaring
            dt = 0.01 * (sigma / vec_len(g)) # 2 tidssteg for å bevege seg 1 sigma
            x_list.append(pos[0])
            y_list.append(pos[1])
            g_list.append(np.sqrt(np.sum(g**2)))
            E_list.append(total_energy(r, t, g))

            F = kin.cpp_kingas.potential_derivative_r(1, r, t) * kin.m0 / (np.prod(kin.mole_weights))
            F_vec = - F * normalize_vec(np.array(pos))
            a = F_vec
            i += 1

            if i > 800 and np.dot(g, pos) < 0:
                break

    return np.array(x_list) / kin.sigma_ij[0, 0], np.array(y_list) / kin.sigma_ij[0, 0], np.array(g_list), np.array(E_list)

def get_path_from_chi(chi, y0=5, xmax=2.5):
    sigma = kin.sigma_ij[0, 0]
    print('chi =', round(chi / np.pi, 2), 'pi')
    r_end = np.sqrt(y0**2 + (xmax + b)**2)
    y_end = - r_end * np.cos(chi)

    x_end = r_end * np.sin(chi) + b
    x = np.array([b, b, x_end]) #/ sigma
    y = np.array([y0, 0, y_end]) #/ sigma
    return x, y

def get_chi_from_path(x, y):
    g_in = np.array([x[1], y[1]]) - np.array([x[0], y[0]])
    g_out = np.array([x[-1], y[-1]]) - np.array([x[-2], y[-2]])
    chi = np.arccos(np.dot(g_in, g_out) / (vec_len(g_in) * vec_len(g_out)))
    if g_out[0] < 0:
        chi = - chi

    print('Chi computed from path :', round(chi / np.pi, 2), 'pi')
    return chi

def get_potential_grid(s_min=0.8, s_max=2.5, Ns=150, Nt=400, ax=None):
    sigma = kin.sigma_ij[0, 0]
    r_list = np.linspace(s_min * sigma, s_max * sigma, Ns)
    t_list = np.linspace(0, 2 * np.pi, Nt)

    x_list = np.linspace(- max(r_list), max(r_list), len(r_list))
    y_list = np.linspace(- max(r_list), max(r_list), len(r_list))

    u_grid = np.full((len(x_list), len(y_list)), np.nan)

    norm = Normalize(vmin=-kin.epsilon_ij[0, 0] / Boltzmann, vmax=kin.epsilon_ij[0, 0] / Boltzmann)

    for r in r_list:
        for t in t_list:
            u = kin.cpp_kingas.potential(1, r, t) / Boltzmann
            x, y = rt_to_xy(r, t)
            x_dist = abs(x_list - x)
            y_dist = abs(y_list - y)

            xi = list(x_dist).index(min(x_dist))
            yi = list(y_dist).index(min(y_dist))

            u_grid[yi, xi] = u

    if ax is None:
        plt.imshow(u_grid, cmap='bwr', norm=norm, extent=[x_list[0], x_list[-1], y_list[0], y_list[-1]])
        plt.show()
    else:
        ax.imshow(u_grid, cmap='bwr', norm=norm, extent=[x_list[0], x_list[-1], y_list[0], y_list[-1]])

def get_force_grid(rmin=0.8,
                   xlim =(-2.5, 2.5),
                   ylim=(-5, 5),
                   N=100, ax=None):
    sigma = kin.sigma_ij[0, 0]

    xmin, xmax = xlim
    ymin, ymax = ylim

    x_list = np.linspace(xmin * sigma, xmax * sigma, N)
    y_list = np.linspace(ymin * sigma, ymax * sigma, N)

    F_grid = np.zeros((len(x_list), len(y_list)))

    for xi, x in enumerate(x_list):
        for yi, y in enumerate(y_list):
            r, t = xy_to_rt(x, y)
            if r < rmin * sigma:
                F = 0
            else:
                F = kin.cpp_kingas.potential_derivative_r(1, r, t)

            F_grid[yi, xi] = F

    F_max = max(F_grid.flatten())
    norm = Normalize(vmin=-F_max, vmax=F_max)

    if ax is None:
        plt.imshow(F_grid, cmap='bwr', norm=norm, extent=[x_list[0]/sigma, x_list[-1]/sigma, y_list[0]/sigma, y_list[-1]/sigma])
        plt.show()
    else:
        ax.imshow(F_grid, cmap='bwr', norm=norm, extent=[x_list[0]/sigma, x_list[-1]/sigma, y_list[0]/sigma, y_list[-1]/sigma])


kin = KineticGas('AR,C1', potential_mode='mie')
sigma = kin.sigma_ij[0, 0]
T, g0, b = 300, 1, 0.3
fig, ax = plt.subplots(1, 1)

get_force_grid(ax=ax, N=100)

x, y, g, _ = get_path(T, g0, b)
g_cmap = get_cmap('plasma')
g_norm = Normalize(vmin=min(g), vmax=max(g))

chi = kin.cpp_kingas.chi(1, T, g0, b * sigma)
x_chi, y_chi = get_path_from_chi(chi)
ax.plot(x_chi, y_chi, color='black', linestyle='--')

print('Minimum distance (numeric):', min(np.sqrt(x**2 + y**2)))
print('Minimum distance (kingas):', kin.cpp_kingas.get_R(1, T, g0, b * sigma) / sigma)
print()
print('Chi (numeric) :', round(get_chi_from_path(x, y) / np.pi, 2), 'pi')
print('Chi (kingas)  :', round(kin.cpp_kingas.chi(1, T, g0, b * sigma) / np.pi, 2), 'pi')

for i in range(1, len(g)):
    ax.plot(x[i - 1:i + 1], y[i - 1:i + 1], color=g_cmap(g_norm(g[i])))
plt.show()


