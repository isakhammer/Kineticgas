from pykingas import py_KineticGas
import os, sys, shutil

args = sys.argv

print('Running __init__.py in', os.path.dirname(__file__), 'with args :')
print(args)
print()
if '-debug' in args or 'debug' in args:
    print('Retrieving debug build')
    shutil.copy(os.path.dirname(__file__)+'/../cpp/debug/KineticGas.so', os.path.dirname(__file__)+'/KineticGas.so')
elif '-release' in args or 'release' in args:
    shutil.copy(os.path.dirname(__file__) + '/../cpp/release/KineticGas.so', os.path.dirname(__file__)+'/KineticGas.so')

KineticGas = py_KineticGas.KineticGas