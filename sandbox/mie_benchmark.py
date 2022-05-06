import numpy as np
import matplotlib.pyplot as plt
from pykingas import KineticGas
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import gas_constant as R
import os

OUT_ROOT = os.path.dirname(__file__) + '/output/'


comps = 'AR,HE'
c1, c2 = comps.split(',')
filename = 'diffusion_test_'+c1+'_'+c2+'_v0'

Vm = 1 / 40
p = 1e5
x1 = 0.5
x = np.array([x1, 1 - x1])

def get_new_version(filename):
    i = 0
    while os.path.isfile(OUT_ROOT + filename + '.csv'):
        i += 1
        filename = ''.join([filename.split('_v')[0], '_v', str(i)])

    return filename

def run(filename):
    filename = get_new_version(filename) + '.csv'
    kin = KineticGas('AR,HE', potential='mie')

    T_list = np.linspace(250, 500)
    D12_list = np.empty_like(T_list)
    for i, T in enumerate(T_list):
        Vm = T * R / p
        print('Computing for T =', T, 'K, p =', p/1e5, 'bar, Vm =', Vm, 'mol/m3')
        D12_list[i] = kin.interdiffusion(T, Vm, x, BH=True, N=1)

    df = pd.DataFrame({'T' : T_list, 'D12' : D12_list})
    df.to_csv(OUT_ROOT + filename)

def plot(filename):
    df = pd.read_csv(OUT_ROOT + filename + '.csv')
    T_data, D12_data = df['T'], df['D12']

    fit_func = lambda T, a, b: a * T**b
    coeff, cov = curve_fit(lambda T, a, b: a * T**b, T_data, D12_data)

    a, b = coeff

    print('Fit is : ', a, 'T^', b, sep='')

    kin_HS = KineticGas(comps, potential='hs')
    D12_hs = np.empty_like(T_data)
    for i, T in enumerate(T_data):
        Vm = T * R / p
        D12_hs[i] = kin_HS.interdiffusion(T, Vm, x, BH=True, N=1)

    plt.scatter(T_data, D12_data, marker='.')
    plt.plot(T_data, fit_func(T_data, a, b))
    plt.plot(T_data, D12_hs)
    plt.show()

filename = get_new_version(filename)
run(filename)
plot(filename)