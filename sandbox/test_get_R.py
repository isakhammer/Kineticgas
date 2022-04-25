from pykingas import KineticGas
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')

kin = KineticGas('AR,C1', potential_mode='mie')
sigma = kin.sigma_ij[0, 0]

func = lambda r: kin.cpp_kingas.get_R_rootfunc(1, T, g, b, r)
derivative = lambda r: kin.cpp_kingas.get_R_rootfunc_derivative(1, T, g, b, r)

T, g, b = 300, 0.1, 1.2 * sigma

r_list = np.linspace(0.99 * sigma, 4 * sigma, 100)
func_list = np.empty_like(r_list)
derivative_list = np.empty_like(r_list)

for i, r in enumerate(r_list):
    func_list[i] = func(r)
    derivative_list[i] = derivative(r)

fig, ax = plt.subplots(1, 1)

r_search_list = []
init_guess_factor = 1
r = init_guess_factor * b
r_search_list.append(r)
r_reset_list = []
next_r = r - func(r) / derivative(r)
r_search_list.append(next_r)
while abs((r - next_r) / sigma) > 1e-5:
    if next_r > 0:
        r = next_r
    else:
        init_guess_factor *= 0.95
        r = init_guess_factor * b
        r_search_list.append(r)
        r_reset_list.append(r)
    next_r = r - func(r) / derivative(r)
    if next_r > 0:
        r_search_list.append(next_r)

print(np.array(r_search_list) / sigma)
ax.plot(np.array(r_search_list) / sigma, [func(r) for r in r_search_list], color='r', marker='o')
plt.plot(np.array(r_reset_list) / sigma, [func(r) for r in r_reset_list], linestyle='', marker='x', color='b')

R = kin.cpp_kingas.get_R(1, T, g, b)

p1, = ax.plot(r_list / sigma, func_list, label='func')
ax.plot(b / sigma, kin.cpp_kingas.get_R_rootfunc(1, T, g, b, b), marker='x', color='b')
ax.plot(R / sigma, func(R), marker='v', color='g')
ax.set_ylabel('func')
ax.set_ylim(-50, 500)
#p2, = twn.plot(r_list / sigma, derivative_list, label='derivative')
#twn.set_ylabel('derivative')

#plt.legend(handles=[p1, p2], loc='center right')
plt.show()