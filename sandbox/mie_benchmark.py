import numpy as np
import matplotlib.pyplot as plt
from pykingas import KineticGas

kin = KineticGas('AR,HE', potential='mie')

T_list = np.linspace(250, 500)
Vm = 1 / 40
x1 = 0.5
x = np.array([x1, 1 - x1])
D12_list = np.empty_like(T_list)
for i, T in enumerate(T_list):
    print(i)
    D12_list[i] = kin.interdiffusion(T, Vm, x, BH=True, N=1)

plt.plot(T_list, D12_list)
plt.show()