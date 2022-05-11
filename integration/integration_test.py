import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')
from scipy.integrate import dblquad
from copy import deepcopy
from pykingas import KineticGas
import time

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

    integ = I.integrate_plane(p1, p2, p3) # V = (1 / 3) * A * h for en pyramide
    if abs(integ - 2) > FLTEPS:
        if do_plot:
            print(1, integ - 2)
            plot_plane(p1, p2, p3)
        return 1

    x = [0, 1, 1]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate_plane(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        if do_plot:
            print(2, integ - 1)
            plot_plane(p1, p2, p3)
        return 2

    x = [1, 1, 2]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate_plane(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        if do_plot:
            plot_plane(p1, p2, p3)
        return 3

    return 0, 0

def expfun(x, y):
    return np.exp(- 0.5 * ((x-5)**2 + (y-5)**2))

def expfun_derx(x, y):
    return - ((x - 3) / 2) * expfun(x, y)

def expfun_dery(x, y):
    return - 2 * (y - 5) * expfun(x, y)

def expfun_dblder(x, y):
    return expfun(x, y) * (((x - 3) / 2)**2 + 4 * (y - 5)**2 - 5 / 2)

def g_integrand(L, r):
    if r == 0:
        return np.exp(- L) * (L**2 + 1)
    else:
        return np.exp(- L) * (L**2) ** (r + 1) + (r + 1) * g_integrand(L, r - 1)

kin = KineticGas('AR,C1', potential='mie')
sigma = kin.sigma_ij[0, 0]
T = 300
N_EVAL = 0
TOT_TIME = 0
def kinfun(g, b, r=2, l=6):
    global N_EVAL, TOT_TIME
    N_EVAL += 1
    t0 = time.process_time()
    val = np.exp(- g ** 2) * g ** (2 * r + 3) * (1 - np.cos(kin.cpp_kingas.chi(1, T, g, b * sigma)) ** l) * b
    t1 = time.process_time()
    TOT_TIME += t1 - t0
    return val


def py_mesh(origin, end, dx, dy, Nxsteps, Nysteps, func, subdomainlimit, f_subdomainlimit, subdomain=False, prev_eval_points=None):
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
        d2fdx2 = (eval_func(x + 2 * dx * abs(Nxsteps), y, Nx + 2 * abs(Nxsteps), Ny) - 2 * eval_func(x + dx * abs(Nxsteps), y, Nx + abs(Nxsteps), Ny) + eval_func(x, y, Nx, Ny)) / (dx * abs(Nxsteps))**2
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

def getpoints(x, y, z):
    return [I.Point(x[i], y[i], z[i]) for i in range(len(x))]

def py_integrate(origin, end, dx, dy, Nxsteps, Nysteps, func, subdomainlimit, f_subdomainlimit, subdomain=False, prev_eval_points=None):
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

    def step(x, y, Nx, Ny, integral):

        x, y = xlist[2], ylist[2]

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
            subintegral, subeval = py_integrate(sub_origin,
                                                  sub_end,
                                                  dx, dy, sub_Nxsteps, sub_Nysteps, func,
                                                  sub_subdomainlimit, sub_f_subdomainlimit,
                                                  subdomain=True,
                                                  prev_eval_points=eval_points)

            for k, v in subeval.items():
                eval_points[k] = v
            integral += subintegral

            x += dx * Nxsteps
            Nx += Nxsteps
            xlist[0] = x
            xlist[1] = x
            xlist[2] = x
            ylist[0] = y
            ylist[1] = y
            ylist[2] = y
            f = eval_func(xlist[2], ylist[2], Nx, Ny)
            zlist[0] = f
            zlist[1] = f
            zlist[2] = f
            return x, y, Nx, Ny, integral

        ylist[0] = ylist[1]
        ylist[1] = ylist[2]
        ylist[2] += dy * Nysteps
        Ny += Nysteps
        zlist[0] = zlist[1]
        zlist[1] = zlist[2]
        zlist[2] = eval_func(xlist[2], ylist[2], Nx, Ny)
        xlist[0] = xlist[1]
        xlist[1] = xlist[2]
        points = getpoints(xlist, ylist, zlist)
        integral += I.integrate_plane(points[0], points[1], points[2])

        xlist[0] = xlist[1]
        xlist[1] = xlist[2]
        xlist[2] += dx * Nxsteps
        Nx += Nxsteps

        ylist[0] = ylist[1]
        ylist[1] = ylist[2]
        ylist[2] -= dy * Nysteps
        Ny -= Nysteps

        zlist[0] = zlist[1]
        zlist[1] = zlist[2]
        zlist[2] = eval_func(xlist[2], ylist[2], Nx, Ny)

        points = getpoints(xlist, ylist, zlist)
        integral += I.integrate_plane(points[0], points[1], points[2])

        return x, y, Nx, Ny, integral

    integral = 0
    x, origin_Nx, y, origin_Ny = origin

    Nx, Ny = origin_Nx, origin_Ny
    Nx_end, Ny_end = end
    xlist = [x, x, x]
    ylist = [y, y, y]
    f = eval_func(x, y, Nx, Ny)
    zlist = [f, f, f]
    x, y, Nx, Ny, integral = step(x, y, Nx, Ny, integral)
    while Ny < Ny_end:

        while origin_Nx < Nx < Nx_end or Nx_end < Nx < origin_Nx:
            x, y, Nx, Ny, integral = step(x, y, Nx, Ny, integral)

        ylist[0] = ylist[1]
        ylist[1] = ylist[2]
        ylist[2] += dy * Nysteps
        Ny += Nysteps
        zlist[0] = zlist[1]
        zlist[1] = zlist[2]
        zlist[2] = eval_func(xlist[2], ylist[2], Nx, Ny)
        xlist[0] = xlist[1]
        xlist[1] = xlist[2]
        points = getpoints(xlist, ylist, zlist)
        integral += I.integrate_plane(points[0], points[1], points[2])
        if Ny < Ny_end:
            Nxsteps *= -1
            x, y, Nx, Ny, integral = step(x, y, Nx, Ny, integral)

    return integral, eval_points

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
    global N_EVAL, TOT_TIME
    if test is True:
        x, y, points = py_mesh((0, 0, 0, 0),
                               (40, 40),
                               0.2, 0.2,
                               4, 4,
                               expfun, 0.6, 0.4)

        fun = expfun

    else:
        N_EVAL, TOT_TIME = 0, 0
        x, y, points = py_mesh((0.1, 0, 0.1, 0), (50, 30), 0.35 / 4, 0.35 / 4, 4, 4, kinfun, 1, 100)
        n_eval = N_EVAL
        tot_time = TOT_TIME
        print('Meshing required', N_EVAL, 'evaluations to generate', len(x), 'mesh points.')
        print('Average evaluation time was :', TOT_TIME / N_EVAL)
        fun = kinfun
    if projection == '3d':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, [fun(xi, yi) for xi, yi in zip(x, y)], marker='.')
        plt.show()
    else:
        plt.plot(x, y, linestyle='', marker='.', color='r')

def test_meshed_integral(projection='2d', test=False):
    global N_EVAL, TOT_TIME
    if test is True:
        integral, points = py_integrate((0, 0, 0, 0),
                           (100, 100),
                           0.1, 0.1,
                           4, 4,
                           expfun, 0.6, 0.4)

        x = []
        y = []
        z = []
        for k, v in points.items():
            x.append(k[0] * 0.2)
            y.append(k[1] * 0.2)
            z.append(v)

    else:
        N_EVAL, TOT_TIME = 0, 0
        integral, points = py_integrate((0.1, 0, 0.1, 0), (50, 30), 0.35 / 4, 0.35 / 4, 4, 4, kinfun, 1, 100)
        n_eval = N_EVAL
        tot_time = TOT_TIME
        x = []
        y = []
        z = []
        for k, v in points.items():
            x.append(k[0] * 0.2)
            y.append(k[1] * 0.2)
            z.append(v)
        print('Meshing required', N_EVAL, 'evaluations to generate', len(x), 'mesh points.')
        print('Average evaluation time was :', TOT_TIME / N_EVAL)
    if projection == '3d':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, marker='.')
        print('Integral is :', integral)
        if test is True:
            a = -np.log(expfun(6, 5))
            print('Exact integral (-inf, inf) is',round(1 / a, 1),'pi. Numeric integral is', integral / np.pi, 'pi')
        plt.show()
    else:
        plt.plot(x, y, linestyle='', marker='.', color='r')

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios' : [1, 0.1]})
#g_grid, b_grid = np.meshgrid(np.linspace(0.3, 5), np.linspace(0.25, 2.5))
#plt.sca(axs[0])


test_meshed_integral(projection='3d', test=False)
plt.show()
#im = plt.pcolormesh(g_grid, b_grid, get_Z(g_grid, b_grid, 2, 6))
#plt.colorbar(im, cax=axs[1])

#plt.contour(xgrid, ygrid, expfun_dblder(xgrid, ygrid))
#plt.sca(axs[0])

#plt.show()

#plot_contour(test=True)
#plot_contour()
#test_plane()
#test_integrate(do_plot=True)
#test_line()
#test_plane()