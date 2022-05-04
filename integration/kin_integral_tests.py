import sys
import numpy as np
from numpy import exp, cos, sqrt
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from cycler import cycler
from itertools import cycle
import mpl_toolkits.mplot3d as a3
plt.style.use('default')
from pykingas import KineticGas
import time
from scipy.integrate import quad

if '-debug' in sys.argv or '-Debug' in sys.argv or '-d' in sys.argv:
    from integration import Integration_d as I
else:
    from integration import Integration_r as I

kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]
T = 300
chi_func = lambda g, b: kin.cpp_kingas.chi(1, T, g, b)

def plot_triangluation(r=1, l=1):
    # Using b = b / sigma
    neval = 0
    def kinfunc(ij, T, g, b, l, r):
        nonlocal neval
        neval += 1
        return 2 * exp(-g**2) * g**(2 * r + 3) * (1 - (cos(chi_func(g, b * sigma)) ** l)) * b
    originxy = (1e-5, 1e-5)
    endxy = (7.5, 5)
    origin = I.Point(originxy[0], originxy[1])
    end = I.Point(endxy[0], endxy[1])
    dx, dy = 0.1, 0.1
    refinement_levels = 8
    subdomain_dblder_limit = 0.05
    verts = I.trisurf(origin, end, dx, dy, refinement_levels, subdomain_dblder_limit, 1, 300, l, r, kinfunc)
    neval = 0
    verts = np.array(verts)
    x, y, z = verts.reshape((verts.shape[0] * verts.shape[1], verts.shape[2])).transpose()
    print('Generated :', verts.shape[0], 'triangles from ', x.shape, 'points')
    cmap = get_cmap('Blues')
    colors = [cmap(i / 3) for i in range(1, 4)]
    cc = cycler(color=colors)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for v, c in zip(verts, cycle(cc)):
        if (abs(v[0, 1] - v[1, 1]) < 1e-10 and abs(v[0, 1] - v[2, 1]) < 1e-10)\
            or (abs(v[0, 0] - v[1, 0]) < 1e-10 and abs(v[0, 0] - v[2, 0]) < 1e-10) :
            continue
        tri = a3.art3d.Poly3DCollection([v])
        tri.set_facecolor(c['color'])
        tri.set_edgecolor('r')
        tri.set(linewidth=0.1)
        ax.add_collection3d(tri)

    ax.set_xlim(originxy[0], endxy[0])
    ax.set_ylim(originxy[1], endxy[1])
    ax.set_zlim(0, 1.25 * max(z))
    ax.set_xlabel(r'$g$ [-]')
    ax.set_ylabel(r'$b$ [$\sigma$]')
    ax.set_title(r'$r$ = '+str(r)+', $\ell$ = '+str(l))
    plt.show()

def integrate_kinfunc(r=1, l=1):
    t0 = time.process_time()
    integral = kin.cpp_kingas.w_spherical(1, T, l, r)
    t = time.process_time() - t0
    print('Integration took', t, 's.')
    print('Cpp-Value is :', integral)
    print('HS integral is :', kin.cpp_kingas.w_HS(1, T, l, r))


def compare_to_HS():
    r_list = [1, 2, 3]
    l_list = [1, 2, 3]
    rlist, llist = np.meshgrid(r_list, l_list)
    numeric = np.empty_like(rlist, float)
    analytic = np.empty_like(rlist, float)
    for ri in range(len(r_list)):
        for li in range(len(l_list)):
            numeric[ri, li] = kin.cpp_kingas.w_spherical(1, T, l_list[li], r_list[ri])
            analytic[ri, li] = kin.cpp_kingas.w_HS(1, T, l_list[li], r_list[ri])

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    print(numeric, analytic)
    ax.scatter(llist, rlist, 100 * (numeric - analytic) / analytic)
    ax.set_xlabel(r'$r$ [-]')
    ax.set_ylabel(r'$\ell$ [-]')
    ax.set_zlabel(r'$\Delta_{HS}W_{r,\ell} / W^{HS}_{r,\ell}$ [%]')
    ax.set_title('Relative deviation between numeric and analytic\ndimentionless collision integrals (%)')
    plt.show()
#integrage_kinfunc(r=1, l=1)
#compare_to_HS()


plot_triangluation(r=1, l=6)
exit(0)
#xs = [1, 1, 3]
#ys = [1, 2, 1]
#zs = [2, 1, 1]
#
#vrts = [xs, ys, zs]
#sx, sy, sz = np.transpose(vrts)
#
#tri = a3.art3d.Poly3DCollection([[xs, ys, zs]])
#tri.set_color('b')
#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')
#ax.add_collection3d(tri)
#ax.scatter(sx, sy, sz, color=['r', 'g', 'b'])
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#plt.show()