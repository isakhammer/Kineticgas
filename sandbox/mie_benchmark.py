import numpy as np
import matplotlib.pyplot as plt
from pykingas import KineticGas
import pandas as pd
from scipy.optimize import curve_fit
from scipy.constants import gas_constant as R
import os, subprocess
from matplotlib.cm import get_cmap
#plt.style.use('default')

OUT_ROOT = os.path.dirname(__file__) + '/output/diffusion/'
PLOT_ROOT = os.path.dirname(__file__) + '/plots/diffusion/'

p = 1e5
x1 = 0.5
x = np.array([x1, 1 - x1])

# For Fuller model
diffusion_volume_map = {'AR' : 16.1, 'HE' : 2.88, 'CO2' : 26.9,
                        'NE' : 5.59, # 'N2' : 17.9, 'O2' : 16.6,
                        'C1' : 16.5 + 4 * 1.98}

mole_weights = np.array([39.95, 4])

mole_weights_map = {'AR' : 39.95, 'HE' : 4, 'CO2' : 44,
                    'NE' : 20.18, # 'N2' : 17.9, 'O2' : 16.6,
                    'C1' : 16}

exp_fit_func = lambda T, a, b: a * T ** b

def get_mole_weights(comps):
    c1, c2 = comps.split(',')
    return np.array([mole_weights_map[c1], mole_weights_map[c2]])

def get_diffusion_volume(comps):
    c1, c2 = comps.split(',')
    return np.array([diffusion_volume_map[c1], diffusion_volume_map[c2]])

def torr_2_Pa(p_torr):
    return p_torr * (1e5/760)

def atm_2_Pa(p_atm):
    return p_atm * 1.01325 * 1e5

def bar_2_atm(p_bar):
    return p_bar / 1.01325

def Pa_2_atm(p_Pa):
    return bar_2_atm(p_Pa / 1e5)

def Fuller(T, p_Pa, comps): #https://pubs.acs.org/doi/pdf/10.1021/ie50677a007?casa_token=i4bNazNv9UAAAAAA:Ie6wz8medLTHGkONA75DQ8N1LHh7utmhcIIhbDTU0LdrLPqaCC81I-hy83BB1UBM1GMi0KH1Dc3ghg
    mA, mB = get_mole_weights(comps)
    vA, vB = get_diffusion_volume(comps)
    return 1e-7 * T**1.75 * np.sqrt((1/mA) + (1/mB)) / (Pa_2_atm(p_Pa) * ((vA**(1/3) + vB**(1/3))**2))

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

    p = 1e5
    T_list = np.linspace(250, 500, 30)
    D12_list = np.empty_like(T_list)
    for i, T in enumerate(T_list):
        Vm = T * R / p
        if i % 5 == 0:
            print('Computing for T =', T, 'K, p =', p/1e5, 'bar, Vm =', Vm, 'mol/m3')
        try:
            D12_list[i] = kin.interdiffusion(T, Vm, x, BH=BH, N=N)
        except Exception as e:
            D12_list[i] = np.nan
            print('Caught exception', e, 'for', kin.comps, ', T, p, N =', T, p, N)

    df = pd.DataFrame({'T' : pd.Series(T_list), 'D12' : pd.Series(D12_list), 'p' : pd.Series(p)})
    df.to_csv(OUT_ROOT + filename)
    print('Saved', OUT_ROOT + filename, '\n')

def run_p(kin, N=1, BH=False):
    filename_root = kin.comps.replace(',', '_') + '_p_' + 'N' + str(N)
    filename = get_new_version(filename_root) # Ensure no overwrite
    print('Running', filename)

    T = 300
    p_list = np.linspace(0.25, 2, 30) * 1e5
    D12_list = np.full_like(p_list, np.nan)
    for i, p in enumerate(p_list):
        Vm = T * R / p
        if i % 5 == 0:
            print('Computing for T =', T, 'K, p =', p/1e5, 'bar, Vm =', Vm, 'mol/m3')
        try:
            D12_list[i] = kin.interdiffusion(T, Vm, x, BH=BH, N=N)
        except Exception as e:
            D12_list[i] = np.nan
            print('Caught exception', e, 'for', kin.comps, ', T, p, N =', T, p, N)

    df = pd.DataFrame({'p' : pd.Series(p_list), 'D12' : pd.Series(D12_list), 'T' : pd.Series(T)})
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
        D12 = df['D12']
        plt.plot(T, D12, marker=markers[i], color=cmap(N / max(N_list)), label=N, markevery=5)

    if any(got_file):
        T = np.ma.MaskedArray(T, mask=np.isnan(D12))
        D12 = np.ma.MaskedArray(D12, mask=np.isnan(D12))
        T = T[~T.mask]
        D12 = D12[~D12.mask]
        coeff, _ = curve_fit(exp_fit_func, T, D12)
        a, b = coeff
        plot_filename = filename_root[:-1] + 'v' + str(v) + '.png'
        D12_fuller = Fuller(T, df['p'][0], comps)
        plt.plot(T, D12_fuller, color='r', label='Fuller')
        plt.ticklabel_format(axis='Y', style='sci', scilimits=(0, 0))
        plt.legend()
        plt.title(comps +'\t' + r'$D_{12} \sim T^{'+str(round(b, 2)) + '}$')
        plt.xlabel(r'$T$ [K]')
        plt.ylabel(r'$D_{12}$ [m$^2$s$^{-1}$]')
        plt.savefig(PLOT_ROOT + plot_filename)
        plt.show()

def plot_p(comps, v=0, N_list=None):
    if N_list is None:
        N_list = [1]
    filename_root = comps.replace(',', '_') + '_p_N'
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

        p = df['p']
        D12 = df['D12']
        plt.plot(p / 1e5, D12, marker=markers[i], color=cmap(N / max(N_list)), label=N, markevery=5)

    if any(got_file):
        p = np.ma.MaskedArray(p, mask=np.isnan(D12))
        D12 = np.ma.MaskedArray(D12, mask=np.isnan(D12))
        p = p[~p.mask]
        D12 = D12[~D12.mask]
        coeff, _ = curve_fit(exp_fit_func, p, D12)
        a, b = coeff
        plot_filename = filename_root[:-1] + 'v' + str(v) + '.png'
        D12_fuller = Fuller(df['T'][0], p, comps)
        plt.plot(p / 1e5, D12_fuller, color='r', label='Fuller')
        plt.ticklabel_format(axis='Y', style='sci', scilimits=(0, 0))
        plt.legend()
        plt.title(comps +'\t' + r'$D_{12} \sim p^{'+str(round(b, 2)) + '}$')
        plt.xlabel(r'$p$ [bar]')
        plt.ylabel(r'$D_{12}$ [m$^2$s$^{-1}$]')
        plt.savefig(PLOT_ROOT + plot_filename)
        plt.show()



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
            plot_p(comps, v=v, N_list=N_list)
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
    run_plotting(v=1, N_list=[1])
    exit(0)
    subprocess.Popen('caffeinate')
    run_computations(BH=False, N_list=[1])
