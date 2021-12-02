#pragma once
#include <vector>
#include "Factorial.h"

class KineticGas{
    public:
    std::vector<double> mole_weights;
    std::vector<std::vector<double>> sigmaij;
    std::vector<double> sigma;
    std::vector<double> mole_fracs;

    const int N; // Degree of approximation
    
    // Binary specific code
    double n, T, m0, M1, M2, x1, x2, m1, m2, sigma1, sigma2, sigma12;

    std::vector<std::vector<double>> A_matrix;
    std::vector<double> delta_vector;

    KineticGas(std::vector<double> init_mole_weights, 
            std::vector<std::vector<double>> init_sigmaij, 
            std::vector<double> init_mole_fracs,
            double init_T,
            double init_p,
            int init_N);

    double omega(int ij, int l, int r);

    std::vector<std::vector<double>> get_A();
    std::vector<double> get_delta();

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