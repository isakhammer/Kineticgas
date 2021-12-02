#pragma once
#include <vector>
#include "Factorial.h"

class KineticGas{
    public:
    std::vector<double> mole_weights;
    std::vector<std::vector<double>> sigmaij;
    std::vector<double> sigma;
    std::vector<double> mole_fracs;
    
    // Binary specific code
    double T, m0, M1, M2, x1, x2, m1, m2, sigma1, sigma2, sigma12;

    KineticGas(std::vector<double> init_mole_weights,
        std::vector<std::vector<double>> init_sigmaij);

    double omega(int ij, int l, int r);

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
double w(int l, int r);

int cpp_tests();

int kingas_tests();

int min(int a, int b);
double min(double a, double b);
int max(int a, int b);
double max(double a, double b);