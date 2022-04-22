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
exit(a)