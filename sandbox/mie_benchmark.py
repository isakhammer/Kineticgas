import numpy as np
import matplotlib.pyplot as plt
from pykingas import KineticGas
import pandas as pd
import os

OUT_ROOT = os.path.dirname(__file__) + '/output/'

comps = 'AR,HE'
c1, c2 = comps.split(',')
filename = 'diffusion_test_'+c1+'_'+c2+'.csv'

kin = KineticGas('AR,HE', potential='mie')

T_list = np.linspace(250, 500)
Vm = 1 / 40
x1 = 0.5
x = np.array([x1, 1 - x1])
D12_list = np.empty_like(T_list)
for i, T in enumerate(T_list):
    print('Computing for T =', T)
    D12_list[i] = kin.interdiffusion(T, Vm, x, BH=True, N=1)

df = pd.DataFrame({'T' : T_list, 'D12' : D12_list})
df.to_csv(OUT_ROOT + filename)

plt.plot(T_list, D12_list)
plt.show()