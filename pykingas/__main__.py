import sys
from pykingas import py_KineticGas

mode = sys.argv[1]
if mode == '-test':
    exit(py_KineticGas.test())
