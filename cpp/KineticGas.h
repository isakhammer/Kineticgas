#pragma once
#include "Factorial.h"
#include "global_params.h"
#include <vector>
#include <map>
#include <functional>


enum potential_modes{
    HS_potential_idx, // Use Hard-sphere potential for omega-integrals
    mie_potential_idx // Use Mie-potential for omega-integrals
};

//typedef  double (KineticGas::*collision_integral)(int ij, int l, int r);  // Pointer to collision integral
class KineticGas{
    public:
    std::vector<double> mole_weights;
    std::vector<std::vector<double>> sigmaij, epsij, la_ij, lr_ij;
    std::vector<double> sigma;
    std::vector<double> mole_fracs;
    const int potential_mode;
    
    double T, m0, M1, M2, x1, x2, m1, m2, 
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

    KineticGas(std::vector<double> init_mole_weights,
        std::vector<std::vector<double>> init_sigmaij,
        std::vector<std::vector<double>> init_epsij,
        std::vector<std::vector<double>> init_la,
        std::vector<std::vector<double>> init_lr,
        int potential_mode);

    // Collision integrals
    double omega(int ij, int l, int r); // Calls the dimentionless collision integral function pointed to by "w_p", selected at initialisation with the "potential_mode" parameter.

    using CollisionIntegralPointer = double(KineticGas::*)(int, double, int, int);
    double w_HS(int ij, double T, int l, int r); // Dimentionless hard-sphere collision integral
    double w_spherical(int ij, double T, int l, int r); // Dimentionless collision integral for spherical potentials
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

    std::vector<std::vector<double>> get_A_matrix(
        double in_T,
        std::vector<double> in_mole_fracs,
        int N);

    std::vector<double> get_delta_vector(
        double T,
        double particle_density,
        int N);
    
    std::vector<std::vector<double>> get_reduced_A_matrix(
        double in_T,
        std::vector<double> in_mole_fracs,
        int N);
    
    std::vector<double> get_alpha_vector(
        double T,
        double particle_density,
        std::vector<double> mole_fracs,
        int N
    );

    double A(int p, int q, int r, int l);
    double A_prime(int p, int q, int r, int l);
    double A_trippleprime(int p, int q, int r, int l);

    double H_ij(int p, int q, int ij);
    double H_i(int p, int q, int ij);
    double H_simple(int p, int q, int i);

    double a(int p, int q);
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