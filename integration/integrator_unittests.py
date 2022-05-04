import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')

if '-debug' in sys.argv or '-Debug' in sys.argv or '-d' in sys.argv:
    from integration import Integration_d as I
else:
    from integration import Integration_r as I

FLTEPS = 1e-12

def test_line(do_plot=False):
    x = [1, 2]
    y = [2, 3]
    p1 = I.Point(x[0], y[0])
    p2 = I.Point(x[1], y[1])
    line = I.get_line(p1, p2)

    r, v = 0, 0
    if abs(line.a - 1) > FLTEPS:
        r, v = 10, line.a
    elif abs(line.b - 1) > FLTEPS:
        r, v = 11, line.b

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        f = lambda x_val: line.a * x_val + line.b
        x_ax = np.linspace(0.9 * x[0], 1.1 * x[1])

        plt.scatter(x, y)
        plt.plot(x_ax, f(x_ax))
        plt.show()

    return r, v

def test_plane(do_plot=False):
    x = [1, 2, 2]
    y = [1, 1, 2]
    z = [0, 6, 0]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    plane = I.get_plane(p1, p2, p3)

    r, v = 0, 0
    if abs(plane.A - 6) > FLTEPS:
        r, v = 20, plane.A
    elif abs(plane.B - (-6)) > FLTEPS:
        r, v = 21, plane.B
    elif abs(plane.C) > FLTEPS:
        r, v = 22, plane.C

    if do_plot is True and (r != 0 or '-force' in sys.argv):
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

    return r, v

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

def test_integrate_plane(do_plot=False):
    x = [0, 1, 2]
    y = [0, 2, 0]
    z = [0, 3, 0]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])

    r, v = 0, 0

    integ = I.integrate_plane(p1, p2, p3) # V = (1 / 3) * A * h for en pyramide
    if abs(integ - 2) > FLTEPS:
        r, v = 30, integ

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        print(30, integ - 2)
        plot_plane(p1, p2, p3)

    x = [0, 1, 1]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate_plane(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        r, v = 31, integ

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        print(31, integ - 1)
        plot_plane(p1, p2, p3)

    x = [1, 1, 2]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate_plane(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        r, v = 32, integ

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        print(32, integ - 1)
        plot_plane(p1, p2, p3)

    return r, v

def test_integrate_2d_linear(do_plot=False):
    integ = I.integrator_test_linear(0, 0, # Origin
                                  10, 10, # End
                                  0.25, 0.25, # dx, dy
                                  4, 0.2)
    r, v = 0, 0
    if abs(integ - 1e3) > FLTEPS: # Should be exact (sans floating point error)
        r, v = 40, integ

    return r, v

def test_integration_expfun(do_plot=False):
    ox, oy = 0, 0
    ex, ey = 10, 10
    dx, dy = 0.2, 0.2
    rlevels, dblderlimit = 4, 0.01
    integ = I.integrator_test(ox, oy, # Origin
                          ex, ey, # End
                          dx, dy, # dx, dy
                          rlevels, dblderlimit) # refinement_levels, dblder_limit

    r, v = 0, 0
    if abs((integ / np.pi) - 1) > 1e-2:
        r, v = 40, integ

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        print('Integration error is :', round(((r / np.pi) - 1) * 100, 2), '%, Integral value is', r / np.pi, 'pi')
        mesh_expfun(ox, oy, # Origin
                    ex, ey, # End
                    dx, dy, # dx, dy
                    rlevels, dblderlimit,
                    projection='3d')

    return r, v

def mesh_expfun(ox, oy, # Origin
                ex, ey, # End
                dx, dy, # dx, dy
                rlevels, dblderlimit, projection='3d'):

    x, y, z = I.mesh_test(ox, oy, # Origin
                          ex, ey, # End
                          dx, dy, # dx, dy
                          rlevels, dblderlimit) # refinement_levels, dblder_limit

    if projection == '3d':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x, y, z, linewidth='0.5')
    else:
        colors = ['r', 'b', 'g', 'black']
        j = 0
        for i in range(1, len(x)):
            plt.plot(x[i-1 : i+1], y[i-1 : i + 1], color=colors[j], linewidth='0.5', marker='o', markersize=j)
            j += 1
            if j == 4:
                j = 0
        #plt.scatter(x, y, marker='o', color='r')
    plt.show()

if __name__ == '__main__':
    do_plot, do_print = False, False
    if '-print' in sys.argv:
        do_print = True
    if '-plot' in sys.argv:
        do_plot = True
    tests = [test_line, test_plane, test_integrate_plane, test_integrate_2d_linear, test_integration_expfun]
    for t in tests:
        r, val = t(do_plot=do_plot)
        if r != 0:
            if do_print:
                print(r, val)
            print('Integration tests failed with exit code :', r)
            exit(r)

    print('Integration tests completed successfully!')
    exit(0)