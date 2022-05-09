/*
    Author : Vegard Gjeldvik Jervell

    Contains :

    Major base-functions for KineticGas. These functions are required for all versions, regardless of what
    potential model they use. Contains summational expressions for the bracket integrals, and iterfaces to
    get the matrices and vectors required to evaluate diffusion coefficients.
*/

#include "KineticGas.h"
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>
#include <math.h>
#include "pybind11/pybind11.h"

#ifdef DEBUG
#define _LIBCPP_DEBUG 1
#endif

#define pprintf(flt) std::printf("%f", flt); std::printf("\n")
#define pprinti(i) std::printf("%i", i); std::printf("\n")

#pragma region // Global helper functions

int min(int a, int b){
    if (a <= b){
        return a;
    }
    return b;
}
double min(double a, double b){
    if (a <= b){
        return a;
    }
    return b;
}
int max(int a, int b){
    if (a >= b){
        return a;
    }
    return b;
}
double max(double a, double b){
    if (a >= b){
        return a;
    }
    return b;
}

int delta(int i, int j){ // Kronecker delta
    if (i == j){
        return 1;
    }
    return 0;
}

std::vector<double> logspace(const double& lmin, const double& lmax, const int& N_gridpoints){
    std::vector<double> grid(N_gridpoints);
    double dx = (lmax - lmin) / (N_gridpoints - 1);
    // Making a logarithmic grid
    // List is inverted when going from lin to log (so that numbers are more closely spaced at the start)
    // Therefore: Count "backwards" to invert the list again so that the smallest r is at r_grid[0]

    // Take a linear grid on (a, b), take log(grid) to get logarithmically spaced grid on (log(a), log(b))
    // Make a linear map from (log(a), log(b)) to (b, a). Because spacing is bigger at the start of the log-grid
    // Map the points on (log(a), log(b)) to (b, a), count backwards such that smaller numbers come first in the returned list
    double A = (lmin - lmax) / log(lmax / lmin);
    double B = lmin - A * log(lmax);
    for (int i = 0; i < N_gridpoints; i++){
        double x = log(lmax - dx * i); // Counting backwards linearly, mapping linear grid to logspace
        grid[i] = A * x + B; // Using linear map from log
    }
    return grid;
}

double erfspace_func(const double& x, const double& lmin, const double& lmax, const double& a, const double& b){
    double r = erf(a * (pow(x, b) - pow(lmin, b)) / (pow(lmax, b) - pow(lmin, b)));
    return r;
}

std::vector<double> erfspace(const double& lmin, const double& lmax, const int& N_gridpoints, double& a, double& b){
    std::vector<double> grid(N_gridpoints);
    double dx = (lmax - lmin) / (N_gridpoints - 1);
    // Making a f(x) grid where f(x) = erf(A * (x^b - lmin) / (lmax^b - lmin))

    double A = (lmin - lmax) / (erfspace_func(lmax, lmin, lmax, a, b) - erfspace_func(lmin, lmin, lmax, a, b));
    double B = lmin - A * erfspace_func(lmax, lmin, lmax, a, b);
    for (int i = 0; i < N_gridpoints; i++){
        double x = lmax - dx * i; // Counting backwards linearly (making linear grid)
        double f = erfspace_func(x, lmin, lmax, a, b); // Mapping linear grid to f-space
        grid[i] = A * f + B; // Linear map from f to lin
    }
    return grid;
}

#pragma endregion

#pragma endregion

#pragma region // Constructor
KineticGas::KineticGas(std::vector<double> init_mole_weights,
        std::vector<std::vector<double>> init_sigmaij,
        std::vector<std::vector<double>> init_epsij,
        std::vector<std::vector<double>> init_la,
        std::vector<std::vector<double>> init_lr,
        int potential_mode)
        : mole_weights{init_mole_weights},
        sigmaij{init_sigmaij},
        epsij{init_epsij},
        la_ij{init_la},
        lr_ij{init_lr},
        m0{0.0},
        potential_mode{potential_mode}
    {

    #ifdef DEBUG
        std::printf("This is a Debug build!\n\n");
    #endif

    for (int i = 0; i < sigmaij.size(); i++){
        sigma.push_back(sigmaij[i][i]);
        m0 += mole_weights[i];
    }
    sigma1 = sigma[0];
    sigma2 = sigma[1];
    sigma12 = sigmaij[0][1];
    sigma_map[1] = sigma1;
    sigma_map[2] = sigma2;
    sigma_map[12] = sigma12;
    sigma_map[21] = sigma12;

    eps1 = epsij[0][0];
    eps2 = epsij[1][1];
    eps12 = epsij[0][1];
    eps_map[1] = eps1;
    eps_map[2] = eps2;
    eps_map[12] = eps12;
    eps_map[21] = eps12;

    la1 = la_ij[0][0];
    la2 = la_ij[1][1];
    la12 = la_ij[0][1];
    la_map[1] = la1;
    la_map[2] = la2;
    la_map[12] = la12;
    la_map[21] = la12;

    lr1 = lr_ij[0][0];
    lr2 = lr_ij[1][1];
    lr12 = lr_ij[0][1];
    lr_map[1] = lr1;
    lr_map[2] = lr2;
    lr_map[12] = lr12;
    lr_map[21] = lr12;

    C1 = (lr1 / (lr1 - la1)) * pow(lr1 / la1, (la1 / (lr1 - la1)));
    C2 = (lr2 / (lr2 - la2)) * pow(lr2 / la2, (la2 / (lr2 - la2)));
    C12 = (lr12 / (lr12 - la12)) * pow(lr12 / la12, (la12 / (lr12 - la12)));
    C_map[1] = C1;
    C_map[2] = C2;
    C_map[12] = C12;
    C_map[21] = C12;

    m1 = mole_weights[0];
    m2 = mole_weights[1];
    M1 = mole_weights[0] / m0;
    M2 = mole_weights[1] / m0;

    w_spherical_integrand_export = std::bind(&KineticGas::w_spherical_integrand, this,
                                            std::placeholders::_1, std::placeholders::_2,
                                            std::placeholders::_3, std::placeholders::_4,
                                            std::placeholders::_5, std::placeholders::_6);

    switch (potential_mode)
    {
    case HS_potential_idx:
        w_p = &KineticGas::w_HS;
        potential_p = &KineticGas::HS_potential;
        p_potential_derivative_r = &KineticGas::HS_potential_derivative;
        p_potential_dblderivative_rr = &KineticGas::HS_potential_dblderivative_rr;
        break;
    case mie_potential_idx:
        w_p = &KineticGas::w_spherical;
        potential_p = &KineticGas::mie_potential;
        p_potential_derivative_r = &KineticGas::mie_potential_derivative;
        p_potential_dblderivative_rr = &KineticGas::mie_potential_dblderivative_rr;
        break;
    default:
        throw "Invalid potential mode!";
    }
}

#pragma endregion

#pragma region // Helper functions

std::vector<std::vector<double>> KineticGas::get_A_matrix(
        const double& T,
        const std::vector<double>& mole_fracs,
        const int& N)
{
    std::vector<std::vector<double>> A_matrix(2 * N + 1, std::vector<double>(2 * N + 1));

    pybind11::gil_scoped_release release;

    std::thread t1(&KineticGas::fill_A_matrix_1, this, std::ref(T), std::ref(mole_fracs), std::ref(N), std::ref(A_matrix));
    std::thread t2(&KineticGas::fill_A_matrix_2, this, std::ref(T), std::ref(mole_fracs), std::ref(N), std::ref(A_matrix));
    std::thread t3(&KineticGas::fill_A_matrix_3, this, std::ref(T), std::ref(mole_fracs), std::ref(N), std::ref(A_matrix));

    for (int p = 1; p <= N; p++){                            // -  -  -  -  -
        for (int q = 0; q <= p; q++){                        // -  -  -  -  -
            A_matrix[p + N][q + N] = a(p, q, T, mole_fracs); // -  -  -  -  -
            A_matrix[q + N][p + N] = A_matrix[p + N][q + N]; // -  -  -  x  -
        }                                                    // -  -  -  x  x
    }

    t1.join(); t2.join(); t3.join();

    pybind11::gil_scoped_acquire acquire;

    return A_matrix;
}

void KineticGas::fill_A_matrix_1( // Fill part of the A-matrix (as indicated by crosses)
        const double& T,                                 // x  -  -  -  -
        const std::vector<double>& mole_fracs,           // x  x  -  -  -
        const int& N,                                    // -  -  -  -  -
        std::vector<std::vector<double>>& A_matrix){     // -  -  -  -  -
                                                         // -  -  -  -  -
    for (int p = - N; p < 0; p++){
        for (int q = - N; q <= - abs(p); q++){
            A_matrix[p + N][q + N] = a(p, q, T, mole_fracs);
            A_matrix[q + N][p + N] = A_matrix[p + N][q + N]; // Matrix is symmetric
        }
    }
}

void KineticGas::fill_A_matrix_2( // Fill part of the A-matrix
        const double& T,                                // -  -  -  -  -
        const std::vector<double>& mole_fracs,          // -  -  -  -  -
        const int& N,                                   // x  x  x  -  -
        std::vector<std::vector<double>>& A_matrix){    // x  x  -  -  -
                                                        // x  -  -  -  -
    for (int p = 0; p <= N; p++){
        for (int q = - N; q <= - abs(p); q++){
            double val = a(p, q, T, mole_fracs);
            A_matrix[p + N][q + N] = val;
            A_matrix[q + N][p + N] = A_matrix[p + N][q + N]; // Matrix is symmetric
        }
    }
}

void KineticGas::fill_A_matrix_3( // Fill part of the A-matrix
        const double& T,                                // -  -  -  -  -
        const std::vector<double>& mole_fracs,          // -  -  -  -  -
        const int& N,                                   // -  -  -  -  -
        std::vector<std::vector<double>>& A_matrix){    // -  -  x  -  -
                                                        // -  x  x  -  -
    for (int p = 1; p <= N; p++){
        for (int q = - p + 1; q < 0; q++){
            A_matrix[p + N][q + N] = a(p, q, T, mole_fracs);
            A_matrix[q + N][p + N] = A_matrix[p + N][q + N]; // Matrix is symmetric
        }
    }
}

std::vector<double> KineticGas::get_delta_vector(
        const double& T,
        const double& particle_density,
        const int& N)
{
    std::vector<double> delta_vector(2 * N + 1);
    delta_vector[N] = (3.0 / (particle_density * 2.0)) * sqrt(2 * BOLTZMANN * T / m0);
    return delta_vector;
}

std::vector<std::vector<double>> KineticGas::get_reduced_A_matrix(
        const double& T,
        const std::vector<double>& mole_fracs,
        const int& N)
{   
    // Get A-matrix, exluding central row and column, where (p == 0 or q == 0)
    std::vector<std::vector<double>> reduced_A(2*N, std::vector<double>(2 * N));
    // Upper left block
    for (int p = - N; p < 0; p++){
        for (int q = - N; q <= p; q++){
            reduced_A[p + N][q + N] = a(p, q, T, mole_fracs);
            reduced_A[q + N][p + N] = reduced_A[p + N][q + N]; // Matrix is symmetric
        }
    }
    //Lower left block (and upper right by symmetry)
    for (int p = 1; p <= N; p++){
        for (int q = - N; q < 0; q++){
            reduced_A[p + N - 1][q + N] = a(p, q, T, mole_fracs);
            reduced_A[q + N][p + N - 1] = reduced_A[p + N - 1][q + N]; // Matrix is symmetric
        }
    }
    //Lower right block
    for (int p = 1; p <= N; p++){
        for (int q = 1; q <= p; q++){
            reduced_A[p + N - 1][q + N - 1] = a(p, q, T, mole_fracs);
            reduced_A[q + N - 1][p + N - 1] = reduced_A[p + N - 1][q + N - 1]; // Matrix is symmetric
        }
    }

    return reduced_A;
}

std::vector<double> KineticGas::get_alpha_vector(
    const double& T,
    const double& particle_density,
    const std::vector<double>& in_mole_fracs,
    const int& N)
{
    std::vector<double> alpha_vector(2 * N);
    alpha_vector[N - 1] = - (15.0 / 4.0) * (in_mole_fracs[1] / particle_density) * sqrt(2 * BOLTZMANN * T / m2);
    alpha_vector[N] = - (15.0 / 4.0) * (in_mole_fracs[0] / particle_density) * sqrt(2 * BOLTZMANN * T / m1);
    return alpha_vector;
}

#pragma endregion

#pragma region // A-functions

double KineticGas::A(const int& p, const int& q, const int& r, const int& l){
    double value{0.0};
    int max_i = min(min(p, q), min(r, p + q + 1 - r));
    for (int i = l - 1; i <= max_i; i++){
        value += ((ipow(8, i) * Fac(p + q - 2 * i) * ipow(-1, l + r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i)) * ipow(4, r)) /
                (Fac(p - i) * Fac(q - i) * Fac(l) * Fac(i + 1 - l) * Fac(r - i) * Fac(p + q + 1 - i - r) * Fac(2 * r + 2)
                    * Fac(p + q + 2 - i) * ipow(4, p + q + 1))) * ((i + 1 - l) * (p + q + 1 - i - r) - l * (r - i));
    }
    return value;
}

double KineticGas::A_prime(const int& p, const int& q, const int& r, const int& l,
                            const double& tmp_M1, const double& tmp_M2){
    double F = (pow(tmp_M1, 2) + pow(tmp_M2, 2)) / (2 * tmp_M1 * tmp_M2);
    double G = (tmp_M1 - tmp_M2) / tmp_M2;

    int max_i = min(p, min(q, min(r, p + q + 1 - r)));
    int max_k;
    int max_w;

    Product p1{1}, p2{1};

    double value{0.0};
    for (int i = l - 1; i <= max_i; i++ ){
        max_w = min(p, min(q, p + q + 1 - r)) - i;
        max_k = min(l, i);
        for (int k = l - 1; k <= max_k; k++){
            for (int w = 0; w <= max_w; w++){
                p1 = ((ipow(8, i) * Fac(p + q - 2 * i - w) * ipow(-1, r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i - w)) * ipow(2, 2 * r) * pow(F, i - k) * pow(G, w)
                        * ((ipow(2, 2 * w - 1) * pow(tmp_M1, i) * pow(tmp_M2, p + q - i - w)) * 2)
                        * (tmp_M1 * (p + q + 1 - i - r - w) * delta(k, l) - tmp_M2 * (r - i) * delta(k, l - 1))
                        ));
                p2 = (Fac(p - i - w) * Fac(q - i - w) * Fac(r - i) * Fac(p + q + 1 - i - r - w) * Fac(2 * r + 2) * Fac(p + q + 2 - i - w) * ipow(4, p + q + 1) * Fac(k) * Fac(i - k) * Fac(w));
                value += p1 / p2; // NB: Gives Bus error if unless v1 and v2 are initialized before adding to value... pls help
            }
        }
    }

    return value;
}

double KineticGas::A_trippleprime(const int& p, const int& q, const int& r, const int& l){
    if (p * q == 0 || l % 2 ){
        return 0.0;
    }
    double value{0.0};
    int max_i = min(p, min(q, min(r, p + q + 1 - r)));
    for (int i = l - 1; i <= max_i; i++){
        value += ((ipow(8, i) * Fac(p + q - (2 * i)) * 2 * ipow(-1, r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i)) * ipow(2, 2 * r) * (((i + 1 - l) * (p + q + 1 - i - r)) - l * (r - i))) /
                    (Fac(p - i) * Fac(q - i) * Fac(l) * Fac(i + 1 - l) * Fac(r - i) * Fac(p + q + 1 - i - r) * Fac(2 * r + 2) * Fac(p + q + 2 - i) * ipow(4, p + q + 1))); 
    }
    value *= pow(0.5, p + q + 1);
    return value;
}
#pragma endregion

#pragma region // H-integrals and a(p, q)
double KineticGas::H_ij(const int& p, const int& q, const int& ij, const double& T){
    double tmp_M1{M1}, tmp_M2{M2};
    if (ij == 21){  // swap indices
        tmp_M1 = M2;
        tmp_M2 = M1;
    }

    double value{0.0};

    int max_l = min(p, q) + 1;
    int max_r;
    for (int l = 1; l <= max_l; l++){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A(p, q, r, l) * omega(12, l, r, T);
        }
    }
    value *= 8 * pow(tmp_M2, p + 0.5) * pow(tmp_M1, q + 0.5);

    return value;
}

double KineticGas::H_i(const int& p, const int& q, const int& ij, const double& T){
    double tmp_M1{M1}, tmp_M2{M2};
    if (ij == 21){  // swap indices
        tmp_M1 = M2;
        tmp_M2 = M1;
    }
    double value{0.0};

    int max_l = min(p, q) + 1;
    int max_r;
    for (int l = 1; l <= max_l; l++){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A_prime(p, q, r, l, tmp_M1, tmp_M2) * omega(12, l, r, T);
        }
    }
    value *= 8;
    
    return value;
}

double KineticGas::H_simple(const int& p, const int& q, const int& i, const double& T){
    double value{0.0};
    int max_l = min(p,q) + 1;
    int max_r;
    for (int l = 2; l <= max_l; l += 2){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A_trippleprime(p, q, r, l) * omega(i, l, r, T);
        }
    }
    value *= 8;
    return value;
}

double KineticGas::a(const int& p, const int& q, const double& T, const std::vector<double>& mole_fracs){
    double x1{mole_fracs[0]}, x2{mole_fracs[1]};
    if (p == 0 || q == 0){
        if (p > 0) return pow(M1, 0.5) * x1 * x2 * H_i(p, q, 12, T);
        else if (p < 0) return - pow(M2, 0.5) * x1 * x2 * H_i(-p, q, 21, T);
        else if (q > 0) return pow(M1, 0.5) * x1 * x2 * H_i(p, q, 12, T);
        else if (q < 0) return - pow(M2, 0.5) * x1 * x2 * H_i(p, -q, 21, T);
        else{  // p == 0 and q == 0
            return M1 * x1 * x2 * H_i(p, q, 12, T);
        }
    }
    else if (p > 0 and q > 0) return pow(x1, 2) * (H_simple(p, q, 1, T)) + x1 * x2 * H_i(p, q, 12, T);

    else if (p > 0 and q < 0) return x1 * x2 * H_ij(p, -q, 12, T);

    else if (p < 0 and q > 0) return x1 * x2 * H_ij(-p, q, 21, T);

    else{  // p < 0 and q < 0
        return pow(x2, 2) * H_simple(-p, -q, 2, T) + x1 * x2 * H_i(-p, -q, 21, T);
    }
}
#pragma endregion

