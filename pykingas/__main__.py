import sys, shutil, os

args = sys.argv

if '-test' in args:
    print('Testing from', __file__)
    try:
        import pyctp
    except:
        print('Missing module dependency pyctp (ThermoPack)')
        
    from pykingas import mie_unittests, KineticGas, py_KineticGas as pykg

    if '-print' in args:
        r = mie_unittests.run_tests(do_print=True)
    elif '-plot' in args:
        r = mie_unittests.run_tests(do_plot=True)
    else:
        r = mie_unittests.run_tests()
    if r != 0:
        print('Mie unittests failed with exit code :', r)
        exit(r)
    elif '-print' not in args and '-plot' not in args:
            print('Mie unittests were successful!')

    kin = KineticGas('AR,HE')
    if '-print' in args:
        r = pykg.test(do_print=True)
    elif '-plot' in args:
        r = pykg.test(plot=True)
    else:
        r = pykg.test()

    if r != 0:
        print('Python test failed with exit code :', r)
    elif '-print' not in args and '-plot' not in args:
        print('Python test was successful!')

    exit(r)
