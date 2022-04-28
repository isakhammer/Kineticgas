import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import dblquad

if '-debug' in sys.argv or '-Debug' in sys.argv or '-d' in sys.argv:
    from integration import Integration_d as I
else:
    from integration import Integration_r as I

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

def test_integrate(do_plot=False):
    x = [0, 1, 1]
    y = [0, 2, 0]
    z = [0, 3, 0]
    p1 = I.Point(x[0], y[0], z[0])
    p2 = I.Point(x[1], y[1], z[1])
    p3 = I.Point(x[2], y[2], z[2])

    plane = I.get_plane(p1, p2, p3)
    integ = I.integrate(p1, p2, p3) # V = (1 / 3) * A * h for en pyramide
    print(integ)
    if do_plot:
        def f(x_val, y_val):
            val = plane.A * x_val + plane.B * y_val + plane.C
            return val

        l12 = I.get_line(p1, p2)
        l13 = I.get_line(p1, p3)
        l23 = I.get_line(p2, p3)
        print(l12.a, l13.a, l23.a)
        y_func = lambda x_val, l: l.a * x_val + l.b

        x_list = np.linspace(0, 2, 10)
        y_list = np.linspace(0, 2, 10)
        x_ax, y_ax = np.meshgrid(x_list, y_list)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, color='r')

        ax.plot_wireframe(x_ax, y_ax, f(x_ax, y_ax), color='black')
        ax.plot(x_list, y_func(x_list, l13), f(x_list, y_func(x_list, l13)), color='r', label='lower')
        ax.plot(x_list, y_func(x_list, l12), f(x_list, y_func(x_list, l12)), color='g', label='upper')
        ax.plot(x_list, y_func(x_list, l23), f(x_list, y_func(x_list, l23)), color='b', label='upper')
        ax.set_ylim(0, 2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.show()


#test_plane()
test_integrate(do_plot=True)
#test_line()
#test_plane()