import sys

args = sys.argv

if '-debug' in args or 'debug' in args:
    from pykingas import KineticGas_d
    __cpp_Module__ = KineticGas_d
else:
    from pykingas import KineticGas_r
    __cpp_Module__ = KineticGas_r

# Expose everything in the __cpp_Module__
for _attr in dir(__cpp_Module__):
    if _attr[:2] != '__': #Exclude macros
        setattr(sys.modules[__name__], _attr, getattr(__cpp_Module__, _attr))

from pykingas import py_KineticGas
KineticGas = py_KineticGas.KineticGas