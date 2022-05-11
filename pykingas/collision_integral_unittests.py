from pykingas import KineticGas, bcolors, suppress_stdout
import matplotlib.pyplot as plt
import numpy as np
import sys

kin = KineticGas('AR,C1', potential='hs')
T = 300

def test_w_vs_HS(do_plot=False, do_print=False):
    # Return value upon failing a test is 200 + 10 * r + l
    # So a failure at r = 3, l = 2 gives exit code 32
    r_list = [1, 2, 3]
    l_list = [1, 2, 3]
    rgrid, lgrid = np.meshgrid(r_list, l_list)
    numeric = np.empty_like(rgrid, float)
    analytic = np.empty_like(rgrid, float)
    print('Computing W-integrals ...')
    r, v = 0, 0
    for ri in range(len(r_list)):
        for li in range(len(l_list)):
            numeric[ri, li] = kin.cpp_kingas.w_spherical(1, T, l_list[li], r_list[ri])
            analytic[ri, li] = kin.cpp_kingas.w_HS(1, T, l_list[li], r_list[ri])

            if abs((numeric[ri, li] / analytic[ri, li]) - 1) > 2.5e-2:
                r = 200 + 10 * rgrid[ri, li] + lgrid[ri, li]
                v = 'For r = ' + str(rgrid[ri, li]) + ', l = ' + str(lgrid[ri, li]) + '\n' \
                    'Numeric HS dimetionless collision integral is : ' + str(numeric[ri, li]) + '\n' \
                    'Analytic HS dimetionless collision integral is : ' + str(analytic[ri, li])
                break
        if r != 0:
            break

    if do_plot is True and (r != 0 or '-force' in sys.argv):
        plot_w_vs_HS(rgrid, lgrid, numeric, analytic)

    return r, v
def plot_w_vs_HS(rgrid, lgrid, numeric, analytic):

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(lgrid, rgrid, 100 * (numeric - analytic) / analytic)
    ax.set_xlabel(r'$r$ [-]')
    ax.set_ylabel(r'$\ell$ [-]')
    ax.set_zlabel(r'$\Delta_{HS}W_{r,\ell} / W^{HS}_{r,\ell}$ [%]')
    ax.set_title('Relative deviation between numeric and analytic\ndimentionless collision integrals (%)')
    plt.show()

def run_tests(do_plot=False, do_print=False):
    '''
        Submodule for testing collision integrals
        Each test in 'tests' must accept two arguments: 'do_plot' and 'do_print' and return two values
        The first is the exit status of the test (0 for successfull, !0 otherwise)
        The second value is some information about the test that failed
    '''
    tests = [test_w_vs_HS]
    if do_plot:
        print('Plotting of mie unittests is not implemented!')
    r = 0
    for t in tests:
        with suppress_stdout('-silent' in sys.argv):
            r, v = t(do_plot, do_print)
        if r != 0:
            if do_print:
                print(r, v)
            print(f'{bcolors.FAIL}Collision integral tests failed with exit code', r, f'{bcolors.ENDC}')
            break
    if r == 0:
        print(f'{bcolors.OKGREEN}Collision integral unittests were successful!{bcolors.ENDC}')
    return r