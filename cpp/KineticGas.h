#pragma once
#include <vector>
#include "Factorial.h"

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
            lr1, lr2, lr12;

    KineticGas(std::vector<double> init_mole_weights,
        std::vector<std::vector<double>> init_sigmaij,
        std::vector<std::vector<double>> init_epsij,
        std::vector<std::vector<double>> init_la,
        std::vector<std::vector<double>> init_lr,
        int potential_mode);

    double omega(int ij, int l, int r);

    using FunctionType = double(KineticGas::*)(int, int, int);
    double omega_HS(int ij, int l, int r);
    double w_HS(int l, int r);
    double omega_Mie(int ij, int l, int r);
    FunctionType omega_p;

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