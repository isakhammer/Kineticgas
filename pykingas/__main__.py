import sys, shutil, os

args = sys.argv

if '-test' in args:
    print('Testing from', __file__)
    try:
        import pyctp
    except:
        print('Missing module dependency pyctp (ThermoPack)')
        
    from pykingas import mie_unittests, KineticGas, py_KineticGas as pykg

    r = mie_unittests.run_tests()
    if r != 0:
        print('Mie unittests failed with exit code :', r)

        if '-print' in args:
            mie_unittests.run_tests(do_print=True)
        if '-plot' in args:
            mie_unittests.run_tests(do_plot=True)

    kin = KineticGas('AR,HE')
    r = pykg.test()
    if r == 0:
        print('Python test was successful!')
    else:
        print('Python test failed with exit code :', r)
    
    if '-print' in args:
        pykg.test(do_print=True)
    if '-plot' in args:
        pykg.test(plot=True)

    exit(r)
