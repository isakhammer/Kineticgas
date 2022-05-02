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

    if do_plot is True:
        f = lambda x_val: line.a * x_val + line.b
        x_ax = np.linspace(0.9 * x[0], 1.1 * x[1])

        plt.scatter(x, y)
        plt.plot(x_ax, f(x_ax))
        plt.show()

    if abs(line.a - 1) > FLTEPS:
        return 10, line.a
    elif abs(line.b - 1) > FLTEPS:
        return 11, line.b

    return 0, 0

def test_plane(do_plot=False):
    x = [1, 2, 2]
    y = [1, 1, 2]
    z = [0, 6, 0]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    plane = I.get_plane(p1, p2, p3)

    if abs(plane.A - 6) > FLTEPS:
        return 20, plane.A
    elif abs(plane.B - (-6)) > FLTEPS:
        return 21, plane.B
    elif abs(plane.C) > FLTEPS:
        return 22, plane.C

    if do_plot is True:
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

    return 0, 0

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

def test_integrate_simple(do_plot=False):
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
        return 30, integ

    x = [0, 1, 1]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate_plane(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        if do_plot:
            print(2, integ - 1)
            plot_plane(p1, p2, p3)
        return 31, integ

    x = [1, 1, 2]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])
    integ = I.integrate_plane(p1, p2, p3)
    if abs(integ - 1) > FLTEPS:
        if do_plot:
            plot_plane(p1, p2, p3)
        return 32, integ

    return 0, 0

def integration_expfun(do_plot=False):
    r = I.test()
    return r

if __name__ == '__main__':
    hold = input('Venter på debugger')
    do_plot, do_print = False, False
    if '-print' in sys.argv:
        do_print = True
    if '-plot' in sys.argv:
        do_plot = True
    tests = [test_line, test_plane, test_integrate_simple, integration_expfun]
    for t in tests:
        r, val = t(do_plot=do_plot)
        if r != 0:
            if do_print:
                print(r, val)
            print('Integration tests failed with exit code :', r)
            exit(r)

    print('Integration tests completed successfully!')
    exit(0)