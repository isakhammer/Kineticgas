import sys, shutil, os

args = sys.argv

if '-test' in args:
    print('Testing from', __file__)
    try:
        import pyctp
    except:
        print('Missing module dependency pyctp (ThermoPack)')
        
    from pykingas import mie_unittests, collision_integral_unittests, KineticGas, py_KineticGas as pykg

    test_pkgs = [mie_unittests.run_tests, collision_integral_unittests.run_tests, pykg.test]

    for test in test_pkgs:
        if '-print' in args and '-plot' in args:
            r = test(do_print=True, do_plot=True)
        elif '-print' in args:
            r = test(do_print=True)
        elif '-plot' in args:
            r = test(do_plot=True)
        else:
            r = test()
        if r != 0:
            exit(r)
    exit(0)
