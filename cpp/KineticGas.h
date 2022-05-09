#pragma once
#include "Factorial.h"
#include "global_params.h"
#include <vector>
#include <map>
#include <functional>
#include <thread>


enum potential_modes{
    HS_potential_idx, // Use Hard-sphere potential for omega-integrals
    mie_potential_idx // Use Mie-potential for omega-integrals
};

// Omega is a function of T, r and l
// Because T is a double, use this struct
struct OmegaPoint{
    int ij, l, r, T_cK;
    OmegaPoint(int ij, int l, int r, double T) : ij{ij}, l{l}, r{r} {
        T_cK = (int) T * 100 + 0.5; // Temperature in cK (10^-2 K)
    };
    bool operator<(const OmegaPoint& other) const {
        if (ij < other.ij) return true;
        else if (ij == other.ij){
            if (l < other.l) return true;
            else if (l == other.l){
                if (r < other.r) return true;
                else if (r == other.r){
                    if (T_cK < other.T_cK) return true;
                }
            }
        }
        return false;
    }

};

//typedef  double (KineticGas::*collision_integral)(int ij, int l, int r);  // Pointer to collision integral
class KineticGas{
    public:
    std::vector<double> mole_weights;
    std::vector<std::vector<double>> sigmaij, epsij, la_ij, lr_ij;
    std::vector<double> sigma;
    const int potential_mode;
    
    double  m0, M1, M2, m1, m2,
            sigma1, sigma2, sigma12, 
            eps1, eps2, eps12,
            la1, la2, la12,
            lr1, lr2, lr12,
            C1, C2, C12;

    std::map<int, double> sigma_map;
    std::map<int, double> eps_map;
    std::map<int, double> la_map;
    std::map<int, double> lr_map;
    std::map<int, double> C_map;
    std::map<OmegaPoint, double> omega_map;

    KineticGas(std::vector<double> init_mole_weights,
        std::vector<std::vector<double>> init_sigmaij,
        std::vector<std::vector<double>> init_epsij,
        std::vector<std::vector<double>> init_la,
        std::vector<std::vector<double>> init_lr,
        int potential_mode);

    // Collision integrals
    double omega(const int& ij, const int& l, const int& r, const double& T); // Calls the dimentionless collision integral function pointed to by "w_p", selected at initialisation with the "potential_mode" parameter.

    using CollisionIntegralPointer = double(KineticGas::*)(const int&, const double&, const int&, const int&);
    double w_HS(const int& ij, const double& T, const int& l, const int& r); // Dimentionless hard-sphere collision integral
    double w_spherical(const int& ij, const double& T, const int& l, const int& r); // Dimentionless collision integral for spherical potentials
    CollisionIntegralPointer w_p; // Will point to one of the above dimentionless collision integrals
    double w_spherical_integrand(const int& ij, const double& T, 
                            const double& g, const double& b, 
                            const int& l, const int& r);
    std::function<double(int, double, double, double, int, int)> w_spherical_integrand_export; // Will bind w_spherical_integrand to this function such that it can be passed to the external integration module
                            
    // Potential models
    double potential(int ij, double r, double theta); // Passes call to the potential corresponding to "potential_mode", using the pointer "potential_p"

    using PotentialPointer = double(KineticGas::*)(int, double, double);
    double HS_potential(int ij, double r, double theta);
    double mie_potential(int ij, double r, double theta);
    PotentialPointer potential_p; // Will point to one of the above potentials

    double potential_derivative_r(int ij, double r, double theta);

    using PotentialDerivativePointer = double(KineticGas::*)(int, double, double);
    double HS_potential_derivative(int ij, double r, double theta);
    double mie_potential_derivative(int ij, double r, double theta);
    PotentialDerivativePointer p_potential_derivative_r; // Will point to one of the above potential derivatives

    double potential_dblderivative_rr(int ij, double r, double theta);

    using PotentialDblDerivativePointer = double(KineticGas::*)(int, double, double);
    double HS_potential_dblderivative_rr(int ij, double r, double theta);
    double mie_potential_dblderivative_rr(int ij, double r, double theta);
    PotentialDblDerivativePointer p_potential_dblderivative_rr; // Will point to one of the above potential derivatives

    // Helper functions for computing dimentionless collision integrals
    double theta(const int ij, const double T, const double g, const double b);
    double theta_lim(const int ij, const double T, const double g);
    double theta_integral(const int ij, const double T, const double R, const double g, const double b);
    double theta_integrand(int ij, double T, double r, double g, double b);
    double theta_integrand_dblderivative(int ij, double T, double r, double g, double b);
    double get_R(int ij, double T, double g, double b);
    double get_R_rootfunc(int ij, double T, double g, double b, double& r);
    double get_R_rootfunc_derivative(int ij, double T, double g, double b, double& r);
    double chi(int ij, double T, double g, double b);
    double chi_HS(int ij, double T, double g, double b);

    std::thread t1;

    void test_th(int& val){
        val += 10;
    }

    std::vector<std::vector<double>> get_A_matrix(
        const double& T,
        const std::vector<double>& in_mole_fracs,
        const int& N);

    void fill_A_matrix_1( // Fill part of the A-matrix
        const double& T,
        const std::vector<double>& mole_fracs,
        const int& N,
        std::vector<std::vector<double>>& A_matrix);

    void fill_A_matrix_2( // Fill another part of the A-matrix
        const double& T,
        const std::vector<double>& mole_fracs,
        const int& N,
        std::vector<std::vector<double>>& A_matrix);

    void fill_A_matrix_3( // Fill another part of the A-matrix
        const double& T,
        const std::vector<double>& mole_fracs,
        const int& N,
        std::vector<std::vector<double>>& A_matrix);

    std::vector<double> get_delta_vector(
        const double& T,
        const double& particle_density,
        const int& N);
    
    std::vector<std::vector<double>> get_reduced_A_matrix(
        const double& T,
        const std::vector<double>& in_mole_fracs,
        const int& N);
    
    std::vector<double> get_alpha_vector(
        const double& T,
        const double& particle_density,
        const std::vector<double>& mole_fracs,
        const int& N
    );

    double A(const int& p, const int& q, const int& r, const int& l);
    double A_prime(const int& p, const int& q, const int& r, const int& l, const double& tmp_M1, const double& tmp_M2);
    double A_trippleprime(const int& p, const int& q, const int& r, const int& l);

    double H_ij(const int& p, const int& q, const int& ij, const double& T);
    double H_i(const int& p, const int& q, const int& ij, const double& T);
    double H_simple(const int& p, const int& q, const int& i, const double& T);

    double a(const int& p, const int& q, const double& T, const std::vector<double>& mole_fracs);
};

int delta(int i, int j);

int cpp_tests();

int kingas_tests();

int min(int a, int b);
double min(double a, double b);
int max(int a, int b);
double max(double a, double b);

double erfspace_func(const double& x, const double& lmin, const double& lmax, const double& a, const double& b);
std::vector<double> erfspace(const double& lmin, const double& lmax, const int& N_gridpoints, double& a, double& b);
std::vector<double> logspace(const double& lmin, const double& lmax, const int& N_gridpoints);