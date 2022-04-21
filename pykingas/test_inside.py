print('Testing')
from pykingas import py_KineticGas as pykg
kin = pykg.KineticGas('AR,HE')
a = pykg.test()
print('Test exit code was :', a)
exit(a)