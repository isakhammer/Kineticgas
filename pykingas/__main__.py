import sys

args = sys.argv[1]
if '-test' in args:
    print('Testing from', __file__)
    try:
        import pyctp
    except:
        print('Missing module dependency pyctp (ThermoPack)')
        
    from pykingas import KineticGas, py_KineticGas as pykg
    kin = KineticGas('AR,HE')
    a = pykg.test()
    if a == 0:
        print('Python test was successful!')
    else:
        print('Python test failed with exit code :', a)
    
    if '-print' in args:
        pykg.test(do_print=True)
    if '-plot' in args:
        pykg.test(plot=True)

    exit(a)
