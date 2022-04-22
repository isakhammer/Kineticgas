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
            lr1, lr2, lr12,
            C1, C2, C12;

    KineticGas(std::vector<double> init_mole_weights,
        std::vector<std::vector<double>> init_sigmaij,
        std::vector<std::vector<double>> init_epsij,
        std::vector<std::vector<double>> init_la,
        std::vector<std::vector<double>> init_lr,
        int potential_mode);

    double omega(int ij, int l, int r);

    using CollisionIntegralPointer = double(KineticGas::*)(int, int, int);
    double w_HS(int ij, int l, int r); // Dimentionless hard-sphere collision integral
    double w_spherical_potential(int ij, int l, int r); // Dimentionless collision integral for spherical potentials
    CollisionIntegralPointer w_p;

    double potential(int ij, double r, double theta);

    using PotentialPointer = double(KineticGas::*)(int, double, double);
    double mie_potential(int ij, double r, double theta);
    double HS_potential(int ij, double r, double theta);
    PotentialPointer potential_p;
    double chi(double g, double b);

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