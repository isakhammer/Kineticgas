from pykingas import KineticGas
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from scipy.constants import Boltzmann
from scipy.integrate import dblquad
plt.style.use('default')

T = 300
kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]



def plot_chi(g=None, b=None):
    b_list = np.linspace(0, 15, 60) * sigma
    g_list = np.linspace(1e-3, 5, 60)

    if g is None and b is None:
        chi_list = np.empty((len(b_list), len(g_list)))
        for bi, b in enumerate(b_list):
            for gi, g in enumerate(g_list):
                chi_list[gi, bi] = kin.cpp_kingas.chi(1, T, g, b)

        b_list, g_list = np.meshgrid(b_list, g_list)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot_wireframe(b_list/sigma, g_list, chi_list / np.pi)
        ax.set_xlabel(r'$b$ [$\sigma$]')
        ax.set_ylabel(r'$g$ [-]')
        ax.set_zlabel(r'$\chi$ [$\pi$]')
        plt.show()
    elif b is None:
        chi_list = np.empty((len(b_list)))
        for bi, b in enumerate(b_list):
            chi_list[bi] = kin.cpp_kingas.theta(1, T, g, b)# * (b / sigma)
        plt.plot(b_list/sigma, chi_list, label=round(g,2))
        plt.plot(b_list / sigma, [chi_list[-1] for _ in b_list], linestyle='--')

    elif g is None:
        chi_list = np.empty((len(g_list)))
        for gi, g in enumerate(g_list):
            chi_list[gi] = kin.cpp_kingas.chi(1, T, g, b)
        plt.plot(g_list, chi_list, label=round(g, 2))
    else:
        print('Either g or b must be None!')

def plot_integrand(l, r, ax, title_l=False, title_r=False):

    b_list = np.linspace(5, 50, 60)
    g_list = np.linspace(1e-5, 5, 60)

    f_list = np.empty((len(b_list), len(g_list)))

    f = lambda g, b: np.exp(- g ** 2) * g ** (2 * r + 3) * (1 - np.cos(kin.cpp_kingas.chi(1, T, g, b * sigma)) ** l) * b

    for bi, b in enumerate(b_list):
        for gi, g in enumerate(g_list):
            f_list[gi, bi] = f(g, b)

    b_list, g_list = np.meshgrid(b_list, g_list)


    ax.plot_wireframe(b_list, g_list, f_list)
    ax.set_xlabel(r'$b$ [$\sigma$]')
    ax.set_ylabel(r'$g$ [-]')
    ax.set_zlabel(r'$f$ [-]')
    if title_r is True:
        ax.set_title('$r$ ='+str(r))
    if title_l is True:
        pass
    #ax.set_title(r'$\ell$ = '+str(l)+', $r$ ='+str(r))

def integrate_quad(l, r): # Bad idea :(
    def integrand(g, b):
        #print(g, b)
        if g < 1e-5 or b * g > 100:
            return 0
        val = np.exp(- g ** 2) * g ** (2 * r + 3) * (1 - np.cos(kin.cpp_kingas.chi(1, T, g, b * sigma)) ** l) * b
        #print(val)
        return val

    I = dblquad(integrand, 0, np.inf, 0, np.inf)
    print(I)

if __name__ == '__main__':
    plot_chi(g=2)
    plt.show()
    exit(0)

    fig = plt.figure()
    for i in range(3):
        for j in range(3):
            print(i, j)
            ax = plt.subplot(3, 3, i*3 + j + 1, projection='3d')
            plot_integrand(i+1, j+1, ax=ax, title_l=(j==0), title_r=(i == 0))

    plt.show()