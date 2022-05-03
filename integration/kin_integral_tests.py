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

kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]
T = 300
chi_func = lambda g, b: kin.cpp_kingas.chi(1, T, g, b)

def mesh_kinfunc(r=2, l=6):
    # Using b = b / sigma
    neval = 0
    def kinfunc(g, b):
        nonlocal neval
        neval += 1
        return exp(-g**2) * g**(2 * r + 3) * (1 - (cos(chi_func(g, b * sigma)) ** l)) * b
    origin = I.Point(1e-5, 1e-5)
    end = I.Point(7.5, 5)
    dx, dy = 0.1, 0.05
    refinement_levels = 8
    subdomain_dblder_limit = 0.05
    t0 = time.process_time()
    mesh = I.mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, kinfunc)
    t = time.process_time() - t0
    x, y, z = mesh
    print('Meshed', len(x), 'points with', neval, 'evaluations in', t, 's.')
    print('Average evaluation time is', t / neval, 's.')
    print('Average time per point is', t / len(x), 's.')
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(x, y, z, marker='.')
    plt.show()

def integrate_kinfunc(r=2, l=6):
    # Using b = b / sigma
    neval = 0
    def kinfunc(g, b):
        nonlocal neval
        neval += 1
        return exp(-g**2) * g**(2 * r + 3) * (1 - (cos(chi_func(g, b * sigma)) ** l)) * b
    origin = I.Point(1e-5, 1e-5)
    end = I.Point(7.5, 5)
    dx, dy = 0.1, 0.05
    refinement_levels = 8
    subdomain_dblder_limit = 0.05
    mesh = I.mesh2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, kinfunc)
    neval = 0
    t0 = time.process_time()
    integral = I.integrate2d(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, kinfunc)
    t = time.process_time() - t0
    x, y, z = mesh
    print('Integrated', len(x), 'points with', neval, 'evaluations in', t, 's.')
    print('Average evaluation time is', t / neval, 's.')
    print('Average time per point is', t / len(x), 's.')
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(x, y, z, marker='.')
    ax.set_title(r'I = '+str(round(integral, 3))+r' $\sigma$')
    plt.show()

integrate_kinfunc()