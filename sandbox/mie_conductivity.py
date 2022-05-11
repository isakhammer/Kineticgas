import numpy as np
import matplotlib.pyplot as plt
from pykingas import KineticGas
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import gas_constant as R
import os, subprocess
from matplotlib.cm import get_cmap
#plt.style.use('default')

OUT_ROOT = os.path.dirname(__file__) + '/output/conductivity/'
PLOT_ROOT = os.path.dirname(__file__) + '/plots/conductivity/'

rho = 40
x1 = 0.5
x = np.array([x1, 1 - x1])

# For Fuller model
diffusion_volume_map = {'CO2' : 16.1, 'C1' : 2.88}

mole_weights = np.array([39.95, 4])

mole_weights_map = {'AR' : 39.95, 'HE' : 4, 'CO2' : 44,
                    'NE' : 20.18, # 'N2' : 17.9, 'O2' : 16.6,
                    'C1' : 16}

def get_new_version(filename_root):
    i = 0
    while os.path.isfile(OUT_ROOT + filename_root + '_v'+str(i) + '.csv'):
        i += 1

    filename = filename_root + '_v'+str(i) + '.csv'

    return filename

def run_T(kin, N=1, BH=False):
    filename_root = kin.comps.replace(',', '_') + '_T_' + 'N' + str(N)
    filename = get_new_version(filename_root) # Ensure no overwrite
    print('Running', filename)

    rho = 40
    T_list = np.linspace(250, 500, 30)
    k_list = np.empty_like(T_list)
    for i, T in enumerate(T_list):
        Vm = 1 / rho
        if i % 5 == 0:
            print('Computing for T =', T, 'K, rho =', rho, 'bar, Vm =', Vm, 'mol/m3')
        try:
            k_list[i] = kin.thermal_conductivity(T, Vm, x, BH=BH, N=N)
        except Exception as e:
            k_list[i] = np.nan
            print('Caught exception', e, 'for', kin.comps, ', T, rho, N =', T, rho, N)

    df = pd.DataFrame({'T' : pd.Series(T_list), 'k' : pd.Series(k_list), 'rho' : pd.Series(rho)})
    df.to_csv(OUT_ROOT + filename)
    print('Saved', OUT_ROOT + filename, '\n')

def run_p(kin, N=1, BH=False):
    filename_root = kin.comps.replace(',', '_') + '_rho_' + 'N' + str(N)
    filename = get_new_version(filename_root) # Ensure no overwrite
    print('Running', filename)

    T = 300
    rho_list = np.linspace(4, 80, 30)
    k_list = np.full_like(rho_list, np.nan)
    for i, rho in enumerate(rho_list):
        Vm = 1 / rho
        if i % 5 == 0:
            print('Computing for T =', T, 'K, rho =', rho, 'bar, Vm =', Vm, 'mol/m3')
        try:
            k_list[i] = kin.thermal_conductivity(T, Vm, x, BH=BH, N=N)
        except Exception as e:
            k_list[i] = np.nan
            print('Caught exception', e, 'for', kin.comps, ', T, rho, N =', T, rho, N)

    df = pd.DataFrame({'rho' : pd.Series(rho_list), 'k' : pd.Series(k_list), 'T' : pd.Series(T)})
    df.to_csv(OUT_ROOT + filename)
    print('Saved', OUT_ROOT + filename, '\n')

def plot_T(comps, v=0, N_list=None):
    if N_list is None:
        N_list = [1]
    filename_root = comps.replace(',', '_') + '_T_N'
    cmap = get_cmap('viridis')
    markers = ['v', 'x', 'o']
    got_file = [True for _ in N_list]
    for i, N in enumerate(N_list):
        filename = filename_root + str(N)+'_v' + str(v)
        data_filename = filename + '.csv'
        try:
            df = pd.read_csv(OUT_ROOT + data_filename)
        except FileNotFoundError:
            got_file[i] = False
            continue

        T = df['T']
        D12 = df['k']
        plt.plot(T, D12, marker=markers[i], color=cmap(N / max(N_list)), label=N, markevery=5)

    if any(got_file):
        plot_filename = filename_root[:-1] + 'v' + str(v) + '.png'
        mod_df = pd.read_csv(OUT_ROOT + 'CO2_C1_T_mod_v0.csv')
        plt.plot(mod_df['T'], mod_df['k'], color='r', label='Model')
        #plt.ticklabel_format(axis='Y', style='sci', scilimits=(0, 0))
        plt.legend()
        plt.title(comps)
        plt.xlabel(r'$T$ [K]')
        plt.ylabel(r'$D_{12}$ [m$^2$s$^{-1}$]')
        #plt.savefig(PLOT_ROOT + plot_filename)
        plt.show()

def plot_rho(comps, v=0, N_list=None):
    if N_list is None:
        N_list = [1]
    filename_root = comps.replace(',', '_') + '_rho_N'
    cmap = get_cmap('viridis')
    markers = ['v', 'x', 'o']
    got_file = [True for _ in N_list]
    for i, N in enumerate(N_list):
        filename = filename_root + str(N)+'_v' + str(v)
        data_filename = filename + '.csv'
        try:
            df = pd.read_csv(OUT_ROOT + data_filename)
        except FileNotFoundError:
            got_file[i] = False
            continue

        rho = df['rho']
        k = df['k']
        plt.plot(rho, k, marker=markers[i], color=cmap(N / max(N_list)), label=N, markevery=5)

    if any(got_file):
        plot_filename = filename_root[:-1] + 'v' + str(v) + '.png'
        mod_df = pd.read_csv(OUT_ROOT + 'CO2_C1_rho_mod_v0.csv')
        plt.plot(mod_df['rho'], mod_df['k'], color='r', label='Model')
        #plt.ticklabel_format(axis='Y', style='sci', scilimits=(0, 0))
        plt.legend()
        plt.title(comps)
        plt.xlabel(r'$\rho$ [mol m$^{-3}$]')
        plt.ylabel(r'$D_{12}$ [m$^2$s$^{-1}$]')
        #plt.savefig(PLOT_ROOT + plot_filename)
        plt.show()


def run_plotting(v=0, N_list=None):
    if N_list is None:
        N_list = [1]
    computed_comps = []
    for k1 in diffusion_volume_map.keys():
        for k2 in diffusion_volume_map.keys():
            comps = k1 + ',' + k2
            if k1 == k2 or comps in computed_comps:
                continue
            plot_T(comps, v=v, N_list=N_list)
            plot_rho(comps, v=v, N_list=N_list)
            computed_comps.append(comps)
            computed_comps.append(k2 + ',' + k1)


def run_computations(BH=False, N_list=None):
    if N_list is None:
        N_list = [1]
    computed_comps = ['HE,CO2', 'CO2,HE']
    for k1 in diffusion_volume_map.keys():
        for k2 in diffusion_volume_map.keys():
            comps = k1 + ',' + k2
            if k1 == k2 or comps in computed_comps:
                continue

            kin = KineticGas(comps, potential='mie')
            for N in N_list:
                print(N, comps)
                run_p(kin, N, BH=BH)
                run_T(kin, N, BH=BH)
            computed_comps.append(comps)
            computed_comps.append(k2 + ',' + k1)

if __name__ == '__main__':
    run_plotting(v=0, N_list=[1])
    exit(0)
    subprocess.Popen('caffeinate')
    run_computations(BH=False, N_list=[1])
    exit(0)

