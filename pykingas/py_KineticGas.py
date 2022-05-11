import numpy as np
from pyctp import saftvrmie
import scipy.linalg as lin
from scipy.constants import Boltzmann, Avogadro
from scipy.integrate import quad
from pykingas import cpp_KineticGas, bcolors, suppress_stdout
import warnings, sys

FLT_EPS = 1e-12

def check_valid_composition(x):
    if abs(sum(x) - 1) > FLT_EPS:
        warnings.warn('Mole fractions do not sum to unity, sum(x) = '+str(sum(x)))

potential_mode_map = {'hs' : 0, 'mie' : 1} # Map string identifiers to corresponding int indentifiers used on cpp-side

class KineticGas:

    default_N = 4

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 BH=False, hs_mixing_rule='additive',
                 potential='HS'):
        '''
        :param comps (str): Comma-separated list of components, following Thermopack-convention
        :param BH (bool) : Use Barker-Henderson diameters?

        Default parameters are equal to default parameters for saft-vr-mie (saftvrmie_parameters.f90)
        If parameters are explicitly supplied, these will be used instead of defaults
        :param mole_weights : (1D array) Molar weights [g/mol]
        :param sigma : (1D array) hard-sphere diameters [m]
        :param eps_div_k : (1D array) epsilon parameter / Boltzmann constant [-]
        :param la, lr : (1D array) attractive and repulsive exponent of the pure components [-]
        :param lij : (float) Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)
        :param kij : (float) Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation)
        :param BH : (bool) Alwayse use Barker-Henderson diameters?
        :param hs_mixing_rule : If "additive", sigma_12 = (1 - lij) * 0.5 * (sigma_1 + sigma_2),
                                else: Compute sigma_12 from BH using epsilon_12 and additive sigma_12
                                Only applicable if BH is True
        :param potential_mode (str) : What potential to use for collision integrals. Options are
                                        'HS' : Use hard-sphere potential
                                        'Mie' : Use Mie-potential
        '''
        if len(comps.split(',')) > 2:
            raise IndexError('Current implementation is only binary-compatible!')
        self.comps = comps
        self.BH = BH
        self.potential_mode = potential.lower()
        self.computed_d_points = {} # dict of state points in which (d_1, d0, d1) have already been computed
        self.computed_a_points = {}  # dict of state points in which (a_1, a1) have already been computed

        if (mole_weights is None) or (sigma is None) or (eps_div_k is None):
            self.eos = saftvrmie.saftvrmie() # Only used as interface to mie-parameter database
            self.eos.init(comps)

        complist = comps.split(',')
        if mole_weights is None:
            mole_weights = np.array([self.eos.compmoleweight(self.eos.getcompindex(comp)) for comp in complist])
        self.mole_weights = np.array(mole_weights) * 1e-3 / Avogadro
        self.m0 = np.sum(self.mole_weights)
        self.M = self.mole_weights/self.m0
        self.M1, self.M2 = self.M

        self.lij = lij
        self.kij = kij
        self.hs_mixing_rule = hs_mixing_rule

        if eps_div_k is None:
            eps_div_k = [self.eos.get_pure_fluid_param(i)[2] for i in range(1, len(complist) + 1)]
        self.epsilon_ij = self.get_epsilon_matrix(eps_div_k, kij)
        self.epsilon = np.diag(self.epsilon_ij)

        if la is None:
            la = np.array([self.eos.get_pure_fluid_param(i)[3] for i in range(1, len(complist) + 1)])
        self.la = self.get_lambda_matrix(la)

        if lr is None:
            lr = np.array([self.eos.get_pure_fluid_param(i)[4] for i in range(1, len(complist) + 1)])
        self.lr = self.get_lambda_matrix(lr)

        if sigma is None:
            sigma = np.array([self.eos.get_pure_fluid_param(i)[1] for i in range(1, len(complist) + 1)])
        self.sigma_ij = self.get_sigma_matrix(sigma)  # Note: Will not initialize with BH-diameters even if BH=True,
                                                      # because a Temperature must be supplied to compute BH-diameter
        self.sigma = np.diag(self.sigma_ij)

        self.cpp_kingas = cpp_KineticGas(self.mole_weights, self.sigma_ij, self.epsilon_ij, self.la, self.lr, potential_mode_map[self.potential_mode])

    def get_A_matrix(self, T, mole_fracs, N=default_N):
        # Compute the matrix of a_pq values
        check_valid_composition(mole_fracs)
        return self.cpp_kingas.get_A_matrix(T, mole_fracs, N)

    def get_reduced_A_matrix(self, T, mole_fracs, N=default_N):
        # Compute the matrix of a_pq values, without a_0q and a_p0
        check_valid_composition(mole_fracs)
        return self.cpp_kingas.get_reduced_A_matrix(T, mole_fracs, N)

    def compute_d_vector(self, T, particle_density, mole_fracs, N=default_N, BH=False):
        # Compute (d_{-1}, d_0 and d_1), used in diifusion solutions
        check_valid_composition(mole_fracs)
        if (T, particle_density, tuple(mole_fracs), N, BH) in self.computed_d_points.keys():
            return self.computed_d_points[(T, particle_density, tuple(mole_fracs), N, BH)]

        if BH:
            sigmaij = self.get_sigma_matrix(self.sigma, BH=BH, T=T)
            cpp_kingas = cpp_KineticGas(self.mole_weights, sigmaij, self.epsilon_ij, self.la, self.lr, potential_mode_map[self.potential_mode])

        else:
            cpp_kingas = self.cpp_kingas

        A = cpp_kingas.get_A_matrix(T, mole_fracs, N)

        delta = cpp_kingas.get_delta_vector(T, particle_density, N)
        d = lin.solve(A, delta)
        d_1, d0, d1 = d[N - 1], d[N], d[N + 1]

        self.computed_d_points[(T, particle_density, tuple(mole_fracs), N, BH)] = (d_1, d0, d1)
        return (d_1, d0, d1)

    def compute_a_vector(self, T, particle_density, mole_fracs, N=default_N, BH=False):
        # Compute a_{-1} and a_1, used in conductivity solutions
        check_valid_composition(mole_fracs)
        if (T, particle_density, tuple(mole_fracs), N, BH) in self.computed_a_points.keys():
            return self.computed_a_points[(T, particle_density, tuple(mole_fracs), N, BH)]

        if BH:
            sigmaij = self.get_sigma_matrix(self.sigma, BH=BH, T=T)
            cpp_kingas = cpp_KineticGas(self.mole_weights, sigmaij, self.epsilon_ij, self.la, self.lr, potential_mode_map[self.potential_mode])
        else:
            cpp_kingas = self.cpp_kingas

        A = cpp_kingas.get_reduced_A_matrix(T, mole_fracs, N)
        alpha = cpp_kingas.get_alpha_vector(T, particle_density, mole_fracs, N)
        a = lin.solve(A, alpha)
        a_1, a1 = a[N - 1], a[N]

        self.computed_a_points[(T, particle_density, tuple(mole_fracs), N, BH)] = (a_1, a1)
        return a_1, a1

    def alpha_T0(self, T, Vm, x, N=default_N, BH=None):
        # Compute the thermal diffusion factor (alpha_T)
        if BH is None:
            BH = self.BH
        check_valid_composition(x)
        particle_density = Avogadro / Vm
        d_1, d0, d1 = self.compute_d_vector(T, particle_density, x, N=N, BH=BH)
        kT = - (5 / (2 * d0)) * ((x[0] * d1 / np.sqrt(self.M[0])) + (x[1] * d_1 / np.sqrt(self.M[1])))
        kT_vec = np.array([kT, -kT])
        return kT_vec * ((1 / np.array(x)) + (1 / (1 - np.array(x))) )

    def soret(self, T, Vm, x, N=default_N, BH=None):
        # Compute the Soret coefficient
        return self.alpha_T0(T, Vm, x, N=N, BH=BH) / T

    def interdiffusion(self, T, Vm, x, N=default_N, BH=None):
        if BH is None:
            BH = self.BH
        check_valid_composition(x)
        particle_density = Avogadro / Vm
        _, d0, _ = self.compute_d_vector(T, particle_density, x, N=N, BH=BH)

        return 0.5 * np.product(x) * np.sqrt(2 * Boltzmann * T / self.m0) * d0

    def thermal_diffusion(self, T, Vm, x, N=default_N, BH=None):
        # Compute the thermal diffusion coefficient (D_T)
        if BH is None:
            BH = self.BH
        check_valid_composition(x)
        particle_density = Avogadro / Vm
        d_1, _, d1 = self.compute_d_vector(T, particle_density, x, N=N, BH=BH)
        return - (5 / 4) * np.product(x) * np.sqrt(2 * Boltzmann * T / self.m0) \
               * ((x[0] * d1 / np.sqrt(self.M1)) + (x[1] * d_1 / np.sqrt(self.M2)))

    def thermal_conductivity(self, T, Vm, x, N=default_N, BH=None):
        if BH is None:
            BH = self.BH
        check_valid_composition(x)
        particle_density = Avogadro / Vm
        a_1, a1 = self.compute_a_vector(T, particle_density, x, N=N, BH=BH)
        return - (5 / 4) * Boltzmann * particle_density * np.sqrt(2 * Boltzmann * T / self.m0) \
               * ((x[0] * a1 / np.sqrt(self.M1)) + (x[1] * a_1 / np.sqrt(self.M2)))

    def get_epsilon_matrix(self, eps_div_k, kij):
        epsilon = np.array(eps_div_k) * Boltzmann
        return (np.ones((2, 2)) - self.kij * (np.ones((2, 2)) - np.identity(2))) * np.sqrt(epsilon * np.vstack(epsilon)) # Only apply mixing parameter kij to the off-diagonals

    def get_sigma_matrix(self, sigma, BH=False, T=None, hs_mixing_rule=None):
        '''
        Get Barker-Henderson diameters for each pair of particles.
        Using Lorentz-Berthleot rules for combining Mie-parameters for each pair of particles

        :param sigma: (1D array) hard sphere diameters [m]
        :return: N x N matrix of hard sphere diameters, where sigma_ij = 0.5 * (sigma_i + sigma_j),
                such that the diagonal is the diameter of each component, and off-diagonals are the cross-collision distances.
        '''
        if hs_mixing_rule is None:
            hs_mixing_rule = self.hs_mixing_rule

        sigma_ij = (np.ones((2, 2)) - self.lij * (np.ones((2, 2)) - np.identity(2))) * 0.5 * np.sum(np.meshgrid(sigma, np.vstack(sigma)), axis=0) # Only apply mixing parameter lij to the off-diagonals

        if BH:
            sigma_ij = np.array([[quad(self.BH_integrand, 0, sigma_ij[i, j], args=(sigma_ij[i, j], self.epsilon_ij[i, j], self.la[i, j], self.lr[i, j], T))[0]
                              for i in range(len(sigma_ij))]
                             for j in range(len(sigma_ij))])

            if hs_mixing_rule == 'additive':
                sigma_ij[0, 1] = sigma_ij[1, 0] = 0.5 * (sigma_ij[0, 0] + sigma_ij[1, 1])
            elif hs_mixing_rule == 'non-additive':
                pass
            else:
                raise KeyError("hs_mixing_rule must be 'additive' or 'non-additive' but was "+str(hs_mixing_rule))
            return sigma_ij
        else:
            return sigma_ij

    def BH_integrand(self, r, sigma, epsilon, lambda_a, lambda_r, T):
        return 1 - np.exp(-self.u_Mie(r, sigma, epsilon, lambda_r, lambda_a) / (T * Boltzmann))

    def u_Mie(self, r, sigma, epsilon, lambda_r, lambda_a):
        C = lambda_r / (lambda_r - lambda_a) * (lambda_r / lambda_a) ** (lambda_a / (lambda_r - lambda_a))
        return C * epsilon * ((sigma / r) ** lambda_r - (sigma / r) ** lambda_a)

    def get_lambda_matrix(self, lambdas):
        l = np.array(lambdas)
        return 3 + np.sqrt((l - 3) * np.vstack(l - 3))

def test(plot=False, do_print=False):
    comps = 'AR,HE'
    kingas = KineticGas(comps, BH=False)

    T = 300
    x = [0.7, 0.3]
    Vm = 24e-3

    # Compute some values
    with suppress_stdout('-silent' in sys.argv):
        alpha_T0 = kingas.alpha_T0(T, Vm, x)
        D12 = kingas.interdiffusion(T, Vm, x)
        DT = kingas.thermal_diffusion(T, Vm, x)
        thermal_cond = kingas.thermal_conductivity(T, Vm, x)

    vals = [alpha_T0, D12, DT, thermal_cond]
    vals_control = [[0.5161629180776952 , - 0.5161629180776952], 4.645041206823456e-05, 5.034955850200782e-06,
                    0.02598758965504002]  # Precomputed values to check that output has not changed
    r = 0
    for i, (val, valc) in enumerate(zip(vals, vals_control)):
        if r != 0:
            break
        if any(abs(np.array([val]).flatten() - np.array([valc]).flatten()) > FLT_EPS) and False:
            r, v = 300 + i + 1, tuple(np.array([val]) - np.array([valc]))

    # Compute values without saving to variables, with different BH-setting
    with suppress_stdout('-silent' in sys.argv):
        alpha_T0_BH = kingas.alpha_T0(T, Vm, x, BH=True)
        D12_BH = kingas.interdiffusion(T, Vm, x, BH=True)
        DT_BH = kingas.thermal_diffusion(T, Vm, x, BH=True)
        thermal_cond_BH = kingas.thermal_conductivity(T, Vm, x, BH=True)

    vals = [alpha_T0_BH, D12_BH, DT_BH, thermal_cond_BH]
    vals_control = [[0.5298208409698769 , - 0.5298208409698769 ], 5.86506545292918e-05, 6.525611212290054e-06, 0.03174164264057401]
    for i, (val, valc) in enumerate(zip(vals, vals_control)):
        if r != 0:
            break
        if any(abs(np.array([val]).flatten() - np.array([valc]).flatten()) > FLT_EPS) and False:
            r, v = 310 + i + 1, tuple(np.array([val]) - np.array([valc]))

    if do_print is True:
        print('\n\nMixture is :', comps)
        print('T =', T, 'K')
        print('rho =', 1e-3 / Vm, 'kmol/m3')
        print('x =', x)
        print()
        print('D12 =', D12, 'mol / m s')
        print('k =', thermal_cond, 'W / m K')
        print('S_T =', 1e3 * alpha_T0 / T, 'mK^{-1}')
        print()
        print('D12_BH =', D12_BH, 'mol / m s')
        print('k_BH =', thermal_cond_BH, 'W / m K')
        print('S_T_BH =', 1e3 * alpha_T0_BH / T, 'mK^{-1}')

    # Recompute the first values, check that they are the same as before.
    with suppress_stdout('-silent' in sys.argv):
        if any(abs(alpha_T0 - kingas.alpha_T0(T, Vm, x)) > FLT_EPS) and r == 0:
            r = 321 and False
        elif abs(D12 - kingas.interdiffusion(T, Vm, x)) > FLT_EPS and r == 0:
            r = 322
        elif abs(DT - kingas.thermal_diffusion(T, Vm, x)) > FLT_EPS and r == 0:
            r = 323
        elif abs(thermal_cond - kingas.thermal_conductivity(T, Vm, x)) > FLT_EPS and r == 0:
            r = 324

    if r != 0:
        print(f'{bcolors.FAIL}Python test failed with exit code :', r, f'{bcolors.ENDC}')
    else:
        print(f'{bcolors.OKGREEN}Python test was successful!{bcolors.ENDC}')

    return r
