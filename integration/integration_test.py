import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')
from scipy.integrate import dblquad
from copy import deepcopy
from pykingas import KineticGas

if '-debug' in sys.argv or '-Debug' in sys.argv or '-d' in sys.argv:
    from integration import Integration_d as I
else:
    from integration import Integration_r as I

FLTEPS = 1e-12

def test_line():
    x = [1, 2]
    y = [2, 3]
    p1 = I.Point(x[0], y[0])
    p2 = I.Point(x[1], y[1])
    line = I.get_line(p1, p2)
    f = lambda x_val: line.a * x_val + line.b

    x_ax = np.linspace(0.9 * x[0], 1.1 * x[1])

    plt.scatter(x, y)
    plt.plot(x_ax, f(x_ax))
    plt.show()

def test_plane():
    x = [1, 2, 2]
    y = [1, 1, 2]
    z = [0, 6, 0]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    plane = I.get_plane(p1, p2, p3)

    f = lambda x_val, y_val: plane.A * x_val + plane.B * y_val + plane.C

    x_ax = np.linspace(0.9 * min(x), 1.1 * max(x))
    y_ax = np.linspace(0.9 * min(y), 1.1 * max(y))
    x_ax, y_ax = np.meshgrid(x_ax, y_ax)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, color='r')

    ax.plot_wireframe(x_ax, y_ax, f(x_ax, y_ax), color='black')
    ax.plot([p1.x, p2.x, p2.x, p2.x, p3.x, p1.x], [p1.y, p2.y, p2.y, p2.y, p3.y, p1.y], [0, 0, p2.z, 0, 0, 0])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()

def plot_plane(p1, p2, p3):
    x = [p1.x, p2.x, p3.x]
    y = [p1.x, p2.y, p3.y]
    z = [p1.z, p2.z, p3.z]

    plane = I.get_plane(p1, p2, p3)

    def f(x_val, y_val):
        val = plane.A * x_val + plane.B * y_val + plane.C
        return val

    l12 = I.get_line(p1, p2)
    l13 = I.get_line(p1, p3)
    l23 = I.get_line(p2, p3)

    y_func = lambda x_val, l: l.a * x_val + l.b

    x_list = np.linspace(min(x), max(x), 10)
    y_list = np.linspace(min(y), max(y), 10)
    x_ax, y_ax = np.meshgrid(x_list, y_list)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, color='r')

    ax.plot_wireframe(x_ax, y_ax, f(x_ax, y_ax), color='black')
    ax.plot(x_list, y_func(x_list, l13), f(x_list, y_func(x_list, l13)), color='r', label='l13')
    ax.plot(x_list, y_func(x_list, l12), f(x_list, y_func(x_list, l12)), color='g', label='l12')
    ax.plot(x_list, y_func(x_list, l23), f(x_list, y_func(x_list, l23)), color='b', label='l23')
    ax.set_ylim(min(y), max(y))
    ax.set_xlim(min(x), max(x))
    ax.set_xlim(min(z), max(z))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()

def test_integrate(do_plot=False):
    x = [0, 1, 2]
    y = [0, 2, 0]
    z = [0, 3, 0]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])

    integ = I.integrate(p1, p2, p3) # V = (1 / 3) * A * h for en pyramide
    if abs(integ - 2) > FLTEPS:
        if do_plot:
            print(1, integ - 2)
            plot_plane(p1, p2, p3)
        return 1

    x = [0, 1, 1]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        if do_plot:
            print(2, integ - 1)
            plot_plane(p1, p2, p3)
        return 2

    x = [1, 1, 2]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        if do_plot:
            plot_plane(p1, p2, p3)
        return 3

    return 0, 0

def expfun(x, y):
    return np.exp(-(((x-3)/2)**2 + (y-5)**2))

def expfun_derx(x, y):
    return - ((x - 3) / 2) * expfun(x, y)

def expfun_dery(x, y):
    return - 2 * (y - 5) * expfun(x, y)

def expfun_dblder(x, y):
    return expfun(x, y) * (((x - 3) / 2)**2 + 4 * (y - 5)**2 - 5 / 2)

kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]
T = 300
N_EVAL = 0
def kinfun(g, b, r=2, l=6):
    global N_EVAL
    N_EVAL += 1
    return np.exp(- g ** 2) * g ** (2 * r + 3) * (1 - np.cos(kin.cpp_kingas.chi(1, T, g, b * sigma)) ** l) * b


def py_mesh(origin, end, dx, dy, Nxsteps, Nysteps, func, subdomainlimit, f_subdomainlimit, subdomain=False, prev_eval_points=None):
    print('Meshing at :', origin[0], origin[2])
    eval_points = {}
    if prev_eval_points is not None:
        for k, v in prev_eval_points.items():
            eval_points[k] = v

    def eval_func(x, y, Nx, Ny):
        if (Nx, Ny) in eval_points.keys():
            return eval_points[(Nx, Ny)]
        else:
            f = func(x, y)
            eval_points[(Nx, Ny)] = f
            return f

    def step(x, y, Nx, Ny):

        d2fdy2 = (eval_func(x, y + 2 * dy * Nysteps, Nx, Ny + 2 * Nysteps) - 2 * eval_func(x, y + dy * Nysteps, Nx, Ny + Nysteps) + eval_func(x, y, Nx, Ny)) / (dy * Nysteps)**2
        d2fdx2 = (eval_func(x + 2 * dx * Nxsteps, y, Nx + 2 * Nxsteps, Ny) - 2 * eval_func(x + dx * Nxsteps, y, Nx + Nxsteps, Ny) + eval_func(x, y, Nx, Ny)) / (dx * Nxsteps)**2
        f = eval_func(x, y, Nx, Ny)
        if (abs(d2fdx2) + abs(d2fdy2) > subdomainlimit or f > f_subdomainlimit) and (abs(Nxsteps) > 1 and abs(Nysteps) > 1):
            sub_origin = (x, Nx, y, Ny)
            sub_Nxsteps = Nxsteps // 2
            sub_Nysteps = Nysteps // 2
            sub_end = (Nx + Nxsteps, Ny + Nysteps)
            sub_subdomainlimit = subdomainlimit * 2
            sub_f_subdomainlimit = f_subdomainlimit * 2
            subx, suby, subeval = py_mesh(sub_origin,
                                          sub_end,
                                          dx, dy, sub_Nxsteps, sub_Nysteps, func,
                                          sub_subdomainlimit, sub_f_subdomainlimit,
                                          subdomain=True,
                                          prev_eval_points=eval_points)

            for k, v in subeval.items():
                eval_points[k] = v
            xlist.extend(subx)
            ylist.extend(suby)

            x += dx * Nxsteps
            Nx += Nxsteps
            xlist.append(x)
            ylist.append(y)
            return x, y, Nx, Ny

        y += dy * Nysteps
        Ny += Nysteps
        xlist.append(x)
        ylist.append(y)

        y -= dy * Nysteps
        Ny -= Nysteps
        x += dx * Nxsteps
        Nx += Nxsteps
        xlist.append(x)
        ylist.append(y)

        return x, y, Nx, Ny

    xlist = []
    ylist = []
    x, origin_Nx, y, origin_Ny = origin

    Nx, Ny = origin_Nx, origin_Ny
    xend, yend = end

    x, y, Nx, Ny = step(x, y, Nx, Ny)
    while Ny < yend:

        while origin_Nx < Nx < xend or xend < Nx < origin_Nx:
            x, y, Nx, Ny = step(x, y, Nx, Ny)

        y += dy * Nysteps
        Ny += Nysteps
        xlist.append(x)
        ylist.append(y)
        if Ny < yend:
            Nxsteps *= -1
            x, y, Nx, Ny = step(x, y, Nx, Ny)

    return deepcopy(xlist), deepcopy(ylist), eval_points

#x, y, points = py_mesh((0, 0, 0, 0), (40, 40), 0.25, 0.25, 4, 4, kinfun)
#print('got points!')
#plt.plot(x, y, linestyle='', marker='.')
def get_Z(g_grid, b_grid, r, l):
    Z = np.empty_like(g_grid)
    for gi, g in enumerate(g_grid):
        for bi, b in enumerate(b_grid):
            #print(gi, g_grid[gi, bi], bi, b_grid[gi, bi])
            Z[gi, bi] = kinfun(g_grid[gi, bi], b_grid[gi, bi], r, l)
    return Z

def plot_cmaps():
    g_grid = np.linspace(0.3, 5, 30)
    b_grid = np.linspace(0.25, 2.5, 30)
    g_grid, b_grid = np.meshgrid(g_grid, b_grid)
    x_grid, y_grid = np.meshgrid(np.linspace(0, 10), np.linspace(0, 10))

    fig, axs = plt.subplots(2, 4, gridspec_kw={'width_ratios':[1, 0.1, 1, 0.1]})

    for i, row in enumerate(axs):
        for j in range(0, len(row), 2):
            r = 2 * i + 2
            l = 2 * j + 2
            Z = get_Z(g_grid, b_grid, r, l)
            maxZ = max(Z.flatten())
            im = row[j].pcolormesh(g_grid, b_grid, Z)
            plt.colorbar(im, cax=row[j + 1])
            #row[j].set_title(round(maxZ, 2))

    axs[1, 0].set_xlabel(r'$g$ [-]')
    axs[1, 2].set_xlabel(r'$g$ [-]')
    axs[0, 0].set_ylabel(r'$b$ [$\sigma$]')
    axs[1, 0].set_ylabel(r'$b$ [$\sigma$]')
    plt.show()

def plot_contour(test=False):
    g_grid = np.linspace(0.3, 5, 30)
    b_grid = np.linspace(0.25, 2.5, 30)
    g_grid, b_grid = np.meshgrid(g_grid, b_grid)

    x_grid, y_grid = np.meshgrid(np.linspace(0, 10), np.linspace(0, 10))

    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios':[1, 0.1]})
    #x, y, points = py_mesh((0.3, 0, 0.25, 0), # Origin
    #                       (50, 40), # Max steps (Nx, Ny)
    #                       0.25 / 4, 0.2 / 4, # dx, dy
    #                       4, 4, # Nxsteps, Nysteps
    #                       expfun, 0.2)



    Z = expfun(x_grid, y_grid)  # get_Z(g_grid, b_grid, 2, 6)
    im = axs[0].contour(x_grid, y_grid, Z)
    plt.colorbar(im, cax=axs[1])

    x, y, points = py_mesh((0, 0, 0, 0),
                           (45, 45),
                           0.2, 0.2,
                           4, 4,
                           expfun, 0.2)

    colors=['r', 'b', 'g', 'black']
    j = 0
    for i in range(len(x)-2):
        axs[0].plot(x[i:i+2], y[i:i+2], linestyle='-', marker='.', color=colors[j], markersize=8)
        j += 1
        if j > 3:
            j = 0
    plt.show()

def plot_mesh(projection='2d', test=False):
    if test is True:
        x, y, points = py_mesh((0, 0, 0, 0),
                               (40, 40),
                               0.2, 0.2,
                               4, 4,
                               expfun, 0.6, 0.4)
        fun = expfun

    else:
        x, y, points = py_mesh((0.3, 0, 0.25, 0), (40, 20), 0.5 / 4, 0.5 / 4, 4, 4, kinfun, 0.2, 100)
        fun = kinfun
    if projection == '3d':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, [fun(xi, yi) for xi, yi in zip(x, y)], marker='.')
        plt.show()
    else:
        plt.plot(x, y, linestyle='', marker='.', color='r')


#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios' : [1, 0.1]})
g_grid, b_grid = np.meshgrid(np.linspace(0.3, 5), np.linspace(0.25, 2.5))
plt.sca(axs[0])
im = plt.pcolormesh(g_grid, b_grid, get_Z(g_grid, b_grid, 2, 6))
plt.colorbar(im, cax=axs[1])
plt.show()
#plt.contour(xgrid, ygrid, expfun_dblder(xgrid, ygrid))

plot_mesh(projection='2d', test=False)

#plot_contour(test=True)
#plot_contour()
#test_plane()
#test_integrate(do_plot=True)
#test_line()
#test_plane()