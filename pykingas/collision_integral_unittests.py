from pykingas import KineticGas
import matplotlib.pyplot as plt
import numpy as np

kin = KineticGas('AR,C1', potential='hs')
T = 300

def test_w_vs_HS(do_plot=False):
    r_list = [1, 2, 3, 4, 5, 6]
    l_list = [1, 2, 3, 4, 5, 6]
    rlist, llist = np.meshgrid(r_list, l_list)
    numeric = np.empty_like(rlist, float)
    analytic = np.empty_like(rlist, float)
    print('Computing W-integrals ...')
    print('-'*36)
    for ri in range(len(r_list)):
        for li in range(len(l_list)):
            print('#', end='')
            numeric[ri, li] = kin.cpp_kingas.w_spherical(1, T, l_list[li], r_list[ri])
            analytic[ri, li] = kin.cpp_kingas.w_HS(1, T, l_list[li], r_list[ri])

            if abs((numeric[ri, li] / analytic[ri, li]) - 1) > 2.5e-2:
                print()
                return 10 * rlist[ri] + l_list[li], (numeric[ri, li], analytic[ri, li])
    print('-'*36)
    if do_plot is True:
        plot_w_vs_HS(rlist, llist, numeric, analytic)

def plot_w_vs_HS(rlist, llist, numeric, analytic):

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(llist, rlist, 100 * (numeric - analytic) / analytic)
    ax.set_xlabel(r'$r$ [-]')
    ax.set_ylabel(r'$\ell$ [-]')
    ax.set_zlabel(r'$\Delta_{HS}W_{r,\ell} / W^{HS}_{r,\ell}$ [%]')
    ax.set_title('Relative deviation between numeric and analytic\ndimentionless collision integrals (%)')
    plt.show()

def run_tests(do_plot=False, do_print=False):
    tests = [test_w_vs_HS]
    if do_plot:
        print('Plotting of mie unittests is not implemented!')
    for t in tests:
        r, val = t(do_plot)
        if r != 0:
            if do_print:
                print(r, val)
                continue
            print('Collision integral tests failed with exit code', r)
            return r
    print('Collision integral unittests were successful!')
    return 0