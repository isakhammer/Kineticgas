import sys
import numpy as np
from numpy import exp, cos, sqrt
import matplotlib.pyplot as plt
plt.style.use('default')
from pykingas import KineticGas
import time

if '-debug' in sys.argv or '-Debug' in sys.argv or '-d' in sys.argv:
    from integration import Integration_d as I
else:
    from integration import Integration_r as I

kin = KineticGas('AR,C1', potential='hs')
sigma = kin.sigma_ij[0, 0]
T = 300
chi_func = lambda g, b: kin.cpp_kingas.chi(1, T, g, b)

def integrate_kinfunc_via_py(r=2, l=6):
    # Using b = b / sigma
    neval = 0
    def kinfunc(ij, T, g, b, l, r):
        nonlocal neval
        neval += 1
        return exp(-g**2) * g**(2 * r + 3) * (1 - (cos(chi_func(g, b * sigma)) ** l)) * b
    origin = I.Point(1e-5, 1e-5)
    end = I.Point(7.5, 5)
    dx, dy = 0.1, 0.05
    refinement_levels = 8
    subdomain_dblder_limit = 0.05
    mesh = I.mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, 1, 300, l, r, kinfunc)
    neval = 0
    t0 = time.process_time()
    integral = I.integrate2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, 1, 300, l, r, kinfunc)
    t = time.process_time() - t0
    x, y, z = mesh
    print('Value is :', integral)
    print('Integrated', len(x), 'points with', neval, 'evaluations in', t, 's.')
    print('Average evaluation time is', t / neval, 's.')
    print('Average time per point is', t / len(x), 's.')
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(x, y, z, linewidth=0.5)
    ax.set_title(r'I = '+str(round(integral, 3))+r' $\sigma$')
    plt.show()

def integrage_kinfunc(r=2, l=6):
    t0 = time.process_time()
    integral = kin.cpp_kingas.w_spherical(1, T, l, r)
    t = time.process_time() - t0
    print('Integration took', t, 's.')
    print('Value is :', integral)
    print('HS integral is :', kin.cpp_kingas.w_HS(1, T, l, r))


integrage_kinfunc()
integrate_kinfunc_via_py()