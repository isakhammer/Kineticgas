#include "KineticGas.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>

#ifdef DEBUG
#define _LIBCPP_DEBUG 1
#endif

constexpr double BOLTZMANN = 1.38064852e-23;
constexpr double GAS_CONSTANT = 8.31446261815324;
constexpr double PI = 3.14159265359;
constexpr double FLTEPS = 1e-10;
#define pprintf(flt) std::printf("%f", flt); std::printf("\n")
#define pprinti(i) std::printf("%i", i); std::printf("\n")

namespace py = pybind11;

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

#pragma region // Tests
int cpp_tests(){
    std::printf("We runnin tests!");
    int r{0};
    r = factorial_tests();
    if (!r) r = kingas_tests();
    return r;
}

int kingas_tests(){
    int d1 = delta(1, 1);
    int d2 = delta(2, 1);
    if (fabs(d1 - 1.0) > FLTEPS || fabs(d2) > FLTEPS){
        return 26;
    }
    
    std::vector<double> Mm{5.5, 10.1};
    std::vector<std::vector<double>> sigmaij {{1.5, 2.0}, {2.0, 2.5}};
    std::vector<std::vector<double>> epsij {{1.5, 2.0}, {2.0, 2.5}};
    std::vector<std::vector<double>> la {{6.0, 6.0}, {6.0, 6.0}};
    std::vector<std::vector<double>> lr {{12.0, 12.0}, {12.0, 12.0}};
    std::vector<double> x {0.3, 0.7};
    KineticGas k{Mm, sigmaij, epsij, la, lr, 0};
    std::vector<double> tsts {k.m0 - 15.6, k.sigma1 - 1.5, k.sigma2 - 2.5, k.sigma12 - 2.0};
    for (double t : tsts){
        if (fabs(t) > FLTEPS){
            return 27;
        }
    }

    double A = k.A(1, 1, 1, 1);
    if (fabs(A) < FLTEPS || fabs(A) > 1e12){
        return 28;
    }
    double A_prime = k.A_prime(1, 1, 1, 1);
    if (fabs(A_prime) < FLTEPS || fabs(A_prime) > 1e12){
        return 29;
    }
    double A_tprime = k.A_trippleprime(2, 3, 2, 2);
    if (fabs(A_tprime) < FLTEPS || fabs(A_tprime) > 1e12){
        return 30;
    }
    return 0;
}
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

    switch (potential_mode)
    {
    case HS_potential_idx:
        w_p = &KineticGas::w_HS;
        potential_p = &KineticGas::HS_potential;
        p_potential_derivative_r = &KineticGas::HS_potential_derivative;
        break;
    case mie_potential_idx:
        w_p = &KineticGas::w_spherical_potential;
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
        double in_T,
        std::vector<double> in_mole_fracs,
        int N)
{
    T = in_T;
    mole_fracs = in_mole_fracs;
    x1 = mole_fracs[0];
    x2 = mole_fracs[1];
    std::vector<std::vector<double>> A_matrix(2*N + 1, std::vector<double>(2 * N + 1));

    for (int p = - N; p <= N; p++){
        for (int q = - N; q <= p; q++){
            A_matrix[p + N][q + N] = a(p, q);
            A_matrix[q + N][p + N] = A_matrix[p + N][q + N]; // Matrix is symmetric
        }
    }

    return A_matrix;
}

std::vector<double> KineticGas::get_delta_vector(
        double in_T,
        double particle_density,
        int N)
{
    std::vector<double> delta_vector(2 * N + 1);
    delta_vector[N] = (3.0 / (particle_density * 2.0)) * sqrt(2 * BOLTZMANN * in_T / m0);
    return delta_vector;
}

std::vector<std::vector<double>> KineticGas::get_reduced_A_matrix(
        double in_T,
        std::vector<double> in_mole_fracs,
        int N)
{   
    // Get A-matrix, exluding central row and column, where (p == 0 or q == 0)
    T = in_T;
    mole_fracs = in_mole_fracs;
    x1 = mole_fracs[0];
    x2 = mole_fracs[1];
    std::vector<std::vector<double>> reduced_A(2*N, std::vector<double>(2 * N));
    // Upper left block
    for (int p = - N; p < 0; p++){
        for (int q = - N; q <= p; q++){
            reduced_A[p + N][q + N] = a(p, q);
            reduced_A[q + N][p + N] = reduced_A[p + N][q + N]; // Matrix is symmetric
        }
    }
    //Lower left block (and upper right by symmetry)
    for (int p = 1; p <= N; p++){
        for (int q = - N; q < 0; q++){
            reduced_A[p + N - 1][q + N] = a(p, q);
            reduced_A[q + N][p + N - 1] = reduced_A[p + N - 1][q + N]; // Matrix is symmetric
        }
    }
    //Lower right block
    for (int p = 1; p <= N; p++){
        for (int q = 1; q <= p; q++){
            reduced_A[p + N - 1][q + N - 1] = a(p, q);
            reduced_A[q + N - 1][p + N - 1] = reduced_A[p + N - 1][q + N - 1]; // Matrix is symmetric
        }
    }

    return reduced_A;
}

std::vector<double> KineticGas::get_alpha_vector(
    double in_T,
    double particle_density,
    std::vector<double> in_mole_fracs,
    int N)
{
    std::vector<double> alpha_vector(2 * N);
    alpha_vector[N - 1] = - (15.0 / 4.0) * (in_mole_fracs[1] / particle_density) * sqrt(2 * BOLTZMANN * in_T / m2);
    alpha_vector[N] = - (15.0 / 4.0) * (in_mole_fracs[0] / particle_density) * sqrt(2 * BOLTZMANN * in_T / m1);
    return alpha_vector;
}

#pragma endregion

#pragma region // A-functions

double KineticGas::A(int p, int q, int r, int l){
    double value{0.0};
    int max_i = min(min(p, q), min(r, p + q + 1 - r));
    for (int i = l - 1; i <= max_i; i++){
        value += ((ipow(8, i) * Fac(p + q - 2 * i) * ipow(-1, l + r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i)) * ipow(4, r)) /
                (Fac(p - i) * Fac(q - i) * Fac(l) * Fac(i + 1 - l) * Fac(r - i) * Fac(p + q + 1 - i - r) * Fac(2 * r + 2)
                    * Fac(p + q + 2 - i) * ipow(4, p + q + 1))) * ((i + 1 - l) * (p + q + 1 - i - r) - l * (r - i));
    }
    return value;
}

double KineticGas::A_prime(int p, int q, int r, int l){
    double F = (pow(M1, 2) + pow(M2, 2)) / (2 * M1 * M2);
    double G = (M1 - M2) / M2;

    int max_i = min(p, min(q, min(r, p + q + 1 - r)));
    int max_k;
    int max_w;

    double value{0.0};
    for (int i = l - 1; i <= max_i; i++ ){
        max_w = min(p, min(q, p + q + 1 - r)) - i;
        max_k = min(l, i);
        for (int k = l - 1; k <= max_k; k++){
            for (int w = 0; w <= max_w; w++){
                value += ((ipow(8, i) * Fac(p + q - 2 * i - w) * ipow(-1, r + i) * Fac(r + 1) * Fac(2 * (p + q + 2 - i - w)) * ipow(2, 2 * r) * pow(F, i - k) * pow(G, w) 
                        * ((ipow(2, 2 * w - 1) * pow(M1, i) * pow(M2, p + q - i - w)) * 2)
                        * (M1 * (p + q + 1 - i - r - w) * delta(k, l) - M2 * (r - i) * delta(k, l - 1))
                        ) / (Fac(p - i - w) * Fac(q - i - w) * Fac(r - i) * Fac(p + q + 1 - i - r - w) * Fac(2 * r + 2) * Fac(p + q + 2 - i - w) * ipow(4, p + q + 1) * Fac(k) * Fac(i - k) * Fac(w))
                        );
            }
        }
    }

    return value;
}

double KineticGas::A_trippleprime(int p, int q, int r, int l){
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
double KineticGas::H_ij(int p, int q, int ij){
    if (ij == 21){  // swap indices
        double tmp{M1};
        M1 = M2;
        M2 = tmp;
    }

    double value{0.0};

    int max_l = min(p, q) + 1;
    int max_r;
    for (int l = 1; l <= max_l; l++){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A(p, q, r, l) * omega(12, l, r);
        }
    }
    value *= 8 * pow(M2, p + 0.5) * pow(M1, q + 0.5);

    
    if (ij == 21){  // swap back
        double tmp{M1};
        M1 = M2;
        M2 = tmp;
    }
    

    return value;
}

double KineticGas::H_i(int p, int q, int ij){
    if (ij == 21){  // swap indices
        double tmp{M1};
        M1 = M2;
        M2 = tmp;
    }
    double value{0.0};

    int max_l = min(p, q) + 1;
    int max_r;
    for (int l = 1; l <= max_l; l++){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A_prime(p, q, r, l) * omega(12, l, r);
        }
    }
    value *= 8;
    
    if (ij == 21){  // swap back
        double tmp{M1};
        M1 = M2;
        M2 = tmp;
    }
    
    return value;
}

double KineticGas::H_simple(int p, int q, int i){
    double value{0.0};
    int max_l = min(p,q) + 1;
    int max_r;
    for (int l = 2; l <= max_l; l += 2){
        max_r = p + q + 2 - l;
        for (int r = l; r <= max_r; r++){
            value += A_trippleprime(p, q, r, l) * omega(i, l, r);
        }
    }
    value *= 8;
    return value;
}

double KineticGas::a(int p, int q){
    if (p == 0 || q == 0){
        if (p > 0) return pow(M1, 0.5) * x1 * x2 * H_i(p, q, 12);
        else if (p < 0) return - pow(M2, 0.5) * x1 * x2 * H_i(-p, q, 21);
        else if (q > 0) return pow(M1, 0.5) * x1 * x2 * H_i(p, q, 12);
        else if (q < 0) return - pow(M2, 0.5) * x1 * x2 * H_i(p, -q, 21);
        else{  // p == 0 and q == 0
            return M1 * x1 * x2 * H_i(p, q, 12);
        }
    }
    else if (p > 0 and q > 0) return pow(x1, 2) * (H_simple(p, q, 1)) + x1 * x2 * H_i(p, q, 12);

    else if (p > 0 and q < 0) return x1 * x2 * H_ij(p, -q, 12);

    else if (p < 0 and q > 0) return x1 * x2 * H_ij(-p, q, 21);

    else{  // p < 0 and q < 0
        return pow(x2, 2) * H_simple(-p, -q, 2) + x1 * x2 * H_i(-p, -q, 21);
    }
}
#pragma endregion

#pragma region // Collision integrals for various potentials

double KineticGas::omega(int ij, int l, int r){
    double w = std::invoke(w_p, this, ij, l, r); // w_p is a pointer to the dimentionless collision integral corresponding to this.potential_mode
    if (ij == 1 || ij == 2){
        return pow(sigma[ij - 1], 2) * sqrt((PI * BOLTZMANN * T) / mole_weights[ij - 1]) * w;
    }
    return 0.5 * pow(sigma12, 2) * sqrt(2 * PI * BOLTZMANN * T / (m0 * M1 * M2)) * w;
}

// Hard-sphere potential
double KineticGas::w_HS(int ij, int l, int r){
    int f = Fac(r + 1).eval();
    if (l % 2 == 0){
        return 0.25 * (2 - ((1.0 / (l + 1)) * 2)) * f;
    }
    return 0.5 * f;
}

// Mie-potential
double KineticGas::w_spherical_potential(int ij, int l, int r){
    throw "Collision integral for Mie potential is not implemented!";
}

#pragma endregion

#pragma region // Various intermolecular potentials

double KineticGas::potential(int ij, double r, double theta){
    return std::invoke(potential_p, this, ij, r, theta);
}

double KineticGas::potential_derivative_r(int ij, double r, double theta){
    return std::invoke(p_potential_derivative_r, this, ij, r, theta);
}

double KineticGas::potential_dblderivative_rr(int ij, double r, double theta){
    return std::invoke(p_potential_dblderivative_rr, this, ij, r, theta);
}

double KineticGas::HS_potential(int ij, double r, double theta){
    if (r > sigma_map[ij]){
        return 0.0;
    }
    return 1e30;
}

double KineticGas::HS_potential_derivative(int ij, double r, double theta){
    if (r < sigma_map[ij]){
        return - 1e30;
    }
    return 0.0;
}

double KineticGas::mie_potential(int ij, double r, double theta){
    return C_map[ij] * eps_map[ij] * (pow(sigma_map[ij] / r, lr_map[ij]) - pow(sigma_map[ij] / r, la_map[ij]));
}

double KineticGas::mie_potential_derivative(int ij, double r, double theta){
    const double lr = lr_map[ij];
    const double la = la_map[ij];
    const double s = sigma_map[ij];
    return C_map[ij] * eps_map[ij] * ((la * pow(s, la) / pow(r, la + 1)) - (lr * pow(s, lr) / pow(r, lr + 1)));
}

double KineticGas::mie_potential_dblderivative_rr(int ij, double r, double theta){
    const double lr = lr_map[ij];
    const double la = la_map[ij];
    const double s = sigma_map[ij];
    return C_map[ij] * eps_map[ij] * ((lr * (lr + 1) * pow(s, lr) / pow(r, lr + 2)) - (la * (la + 1) * pow(s, la) / pow(r, la + 2)));
}

#pragma endregion

#pragma region // Helper funcions for computing dimentionless collision integrals

double KineticGas::theta(const int ij, const double T, const double r_prime, const double g, const double b, const int init_N_gridpoints){
    double lower_limit = r_prime;
    double upper_limit = lower_limit;


    double eps_tol = 1e-3; // Error tolerance
    double convergence_threshold = 0.99; // Integral after 90% of steps must be >= convergence_threshold * integral at 100% of steps
    double eps_0; // Upper limit for error in the first integration step (not tolerance, but mathematical maximum error)
    double eps_N; // Upper limit for error in the last integratio step (not tolerance, but mathematical maximum error)
    double rel_eps_tot = 0; // Upper limit for total relative error (not tolerance, but mathematical maximum error)

    double total_integral = 0;
    double integral_09; // Integral at 90% completion
    double A_coeff, B_coeff; // Coefficients for y = Ax + B, that interpolates the integrand in points (r1, r2)

    int N_integrated_points = 0; // To keep track of total number of integrated points, accross several iterations of the outermost loop
    double upper_limit_cutoff_factor = 1e-5;
    double theta_0 = theta_integrand(ij, T, lower_limit, g, b);

    do{ // Compute integral, recompute if not converged
        // Until convergence (by adjusting upper_limit_cutoff_factor)
        std::vector<double> r_grid;

        lower_limit = upper_limit; // On iterations > 1, this prevents the integral up to "upper_limit" from being recomputed.

        do {
            upper_limit += lower_limit;
        } while (theta_integrand(ij, T, upper_limit, g, b) > upper_limit_cutoff_factor * theta_0);

        int N_gridpoints = 0;

        do{ // Increase number of integration points
            // Until Error in the last integratio step is below tolerance (By adjusting N_gridpoints)

            N_gridpoints += init_N_gridpoints;
            double erfspace_a = 3; // Will increment to 4.5 immediately
            double erfspace_b = 0.9; // Will increment to 0.81 immediately

            do { // Increase density of points at start of integral
                 // Until Error in the first integration step is below tolerance (by adjusting erfspace_a and erfspace_b)
                erfspace_a *= 1.5;
                erfspace_b *= 0.9;
                r_grid = erfspace(lower_limit, upper_limit, N_gridpoints, erfspace_a, erfspace_b);

                double t0 = theta_integrand(ij, T, r_grid[0], g, b);
                double t1 = theta_integrand(ij, T, r_grid[1], g, b);
                double d2tdr2 = theta_integrand_dblderivative(ij, T, r_grid[0], g, b); // Second derivative is strictly decreasing, so max is at first integration point

                double integral = 0;
                A_coeff = (t1 - 0) / (r_grid[0] - r_grid[1]);
                B_coeff = t0 - A_coeff * r_grid[0];
                integral += A_coeff * (pow(r_grid[0], 2) - pow(r_grid[1], 2)) / 2 + B_coeff * (r_grid[0] - r_grid[1]);

                eps_0 = d2tdr2 * pow(r_grid[1] - r_grid[0], 3) / 12.0; // Upper limit for error

                #ifdef DEBUG
                    if (eps_0 > eps_tol * N_gridpoints || isnan(eps_0)){
                            std::printf("Error in first step is less than : %E\n", eps_0 / N_gridpoints);
                            std::printf("Adjusting A, B to : %E, %E\n\n", erfspace_a * 1.5, erfspace_b * 0.9);
                    }
                #endif

            } while (eps_0 > eps_tol * N_gridpoints);

            const int N = N_gridpoints - 1;
            double tN_1 = theta_integrand(ij, T, r_grid[N - 1], g, b);
            double tN = theta_integrand(ij, T, r_grid[N], g, b);
            double d2tdr2 = theta_integrand_dblderivative(ij, T, r_grid[N - 1], g, b); // Second derivative is strictly decreasing, so max is at first integration point

            double integral = 0;
            A_coeff = (tN - tN_1) / (r_grid[N] - r_grid[N - 1]);
            B_coeff = tN - A_coeff * r_grid[N];
            integral += A_coeff * (pow(r_grid[N], 2) - pow(r_grid[N - 1], 2)) / 2 + B_coeff * (r_grid[N] - r_grid[N - 1]);

            eps_N = d2tdr2 * pow(r_grid[N] - r_grid[N - 1], 3) / 12;

            #ifdef DEBUG
                if (eps_N > eps_tol * N_gridpoints || isnan(eps_N)){
                        std::printf("Error in last step is less than : %E\n", eps_N / (N_gridpoints - 50));
                        std::printf("Error in first step is less than : %E\n", eps_0 / (N_gridpoints - 50));
                        std::printf("Adjusting N to : %i\n\n", N_gridpoints);
                }
            #endif

        } while (eps_N > eps_tol * N_gridpoints);

        // Trapezoid rule (piecewise linear interpolation) integration
        double integrand1;
        double integrand2;
        double r1 = r_grid[0];
        double r2 = r_grid[1];
        integrand1 = theta_integrand(ij, T, r1, g, b);
        integrand2 = theta_integrand(ij, T, r2, g, b);
        A_coeff = (integrand2 - integrand1) / (r2 - r1);
        B_coeff = integrand1 - A_coeff * r1;
        total_integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1);
        rel_eps_tot += theta_integrand_dblderivative(ij, T, r1, g, b) * pow(r2 - r1, 3) / 12.0;
        // N_integrated_points is the number of points from previous iterations
        int i; // Loop is split in two parts, so "i" must be initialized outside loop
        std::printf("Printing (i, r [sigma], theta, eps\n");
        for (i = 2; (N_integrated_points + i) < 9 * (N_integrated_points + N_gridpoints) / 10; i++){ // Integrate up to 90% of the total integration steps
            r1 = r2;
            integrand1 = integrand2;
            r2 = r_grid[i];
            integrand2 = theta_integrand(ij, T, r2, g, b);
            A_coeff = (integrand2 - integrand1) / (r2 - r1);
            B_coeff = integrand1 - A_coeff * r1;
            total_integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1);
            rel_eps_tot += theta_integrand_dblderivative(ij, T, r1, g, b) * pow(r2 - r1, 3) / 12.0;
            std::printf("%i \t %E \t %E \t %E\n", i + N_integrated_points, r1 / sigma_map[ij], total_integral, rel_eps_tot);
        }
        integral_09 = total_integral; // remember value at 90% to check convergence after finishing
        for (;i < N_gridpoints; i++){ // Finish integral
            r1 = r_grid[i - 1];
            r2 = r_grid[i];
            integrand1 = integrand2;
            integrand2 = theta_integrand(ij, T, r2, g, b);
            A_coeff = (integrand2 - integrand1) / (r2 - r1);
            B_coeff = integrand1 - A_coeff * r1;
            total_integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1);
            rel_eps_tot += theta_integrand_dblderivative(ij, T, r1, g, b) * pow(r2 - r1, 3) / 12.0;
            std::printf("%i \t %E \t %E \t %E\n", i + N_integrated_points, r1 / sigma_map[ij], total_integral, rel_eps_tot);
        }
        std::printf("\n\n");
        N_integrated_points += N_gridpoints;

        upper_limit_cutoff_factor *= 0.1;
        #ifdef DEBUG
            if (integral_09 < convergence_threshold * total_integral || isnan(total_integral)){
                std::printf("theta = %E pi\n", total_integral / PI);
                std::printf("Change in final 10%% is : %E %% of final value\n", (1 - integral_09 / total_integral) * 100.0);
                std::printf("Adjusting cutoff factor to : %E\n", upper_limit_cutoff_factor);
                std::printf("N_gridpoints and r_grid.size() are %i, %i\n\n", N_gridpoints, r_grid.size());
            }
        #endif

    } while (integral_09 < convergence_threshold * total_integral);

    rel_eps_tot = rel_eps_tot / total_integral;

    #ifdef DEBUG
        std::printf("For b = %E sigma, g = %E\n", b / sigma_map[ij], g);
        std::printf("init, start, end (sigma) = %E %E, %E\n", r_prime / sigma_map[ij], lower_limit / sigma_map[ij], upper_limit / sigma_map[ij]);
        std::printf("Total gridpoints = %i\n", N_integrated_points);
        std::printf("Computed theta = %E pi\n", total_integral / PI);
        std::printf("Relative error is less than %E %% \n\n", rel_eps_tot * 100);
    #endif

    return total_integral;
}

double KineticGas::theta_integrand(int ij, double T, double r, double g, double b){
    // Passing dummy value "1.0" to potential. Mie potential is spherical, so not a function of last parameter (theta).
    return pow((pow(r, 4) / pow(b, 2)) * (1.0 - potential(ij, r, 1.0) / (BOLTZMANN * T * pow(g, 2))) - pow(r, 2), -0.5);
}

double KineticGas::theta_integrand_dblderivative(int ij, double T, double r, double g, double b){
    // Expressing the integrand as f = (core)^{-1/2}
    const double a = 1.0 / (pow(b, 2) * BOLTZMANN * T * pow(g, 2));
    const double u = potential(ij, r, 1.0);
    const double u_prime = potential_derivative_r(ij, r, 1.0);
    const double u_dblprime = potential_dblderivative_rr(ij, r, 1.0);
    const double core = pow(r, 4) / pow(b, 2) - a * pow(r, 4) * u - pow(r, 2);
    const double core_prime = 4 * pow(r, 3) / pow(b, 2) - a * (4 * pow(r, 3) * u + pow(r, 4) * u_prime) - 2 * r;
    const double core_dblprime = 12.0 * pow(r, 2) / pow(b, 2) - a * (12 * pow(r, 2) * u + 8 * pow(r, 3) * u_prime + pow(r, 4) * u_dblprime) - 2;

    return (3.0 / 4.0) * pow(core, -2.5) * core_prime - 0.5 * pow(core, - 1.5) * core_dblprime;
}

double KineticGas::get_R_rootfunc(int ij, double T, double g, double b, double& r){
    return (potential(ij, r, 1.0) / (BOLTZMANN * T * pow(g, 2))) + pow(b / r, 2) - 1;
}

double KineticGas::get_R_rootfunc_derivative(int ij, double T, double g, double b, double& r){
    return (potential_derivative_r(ij, r, 1.0) / (BOLTZMANN * T * pow(g, 2))) - 2 * pow(b, 2) / pow(r, 3);
}

double KineticGas::get_R(int ij, double T, double g, double b){
    // Newtons method
    double tol = 1e-5; // Relative to sigma_map[ij]
    double init_guess_factor = 1.0;
    double r = init_guess_factor * b;
    double f = get_R_rootfunc(ij, T, g, b, r);
    double dfdr = get_R_rootfunc_derivative(ij, T, g, b, r);
    double next_r = r - f / dfdr;
    while (abs((r - next_r) / sigma_map[ij]) > tol){
        if (next_r < 0){
            init_guess_factor *= 0.95;
            r = init_guess_factor * b;
            #ifdef DEBUG
                std::printf("Initial guess for R failed (r < 0), reducing to %E sigma\n\n", r / sigma_map[ij]);
            #endif
        }
        else if (f < 0 && f / dfdr < 0){
            init_guess_factor *= 0.95;
            r = init_guess_factor * b;
            #ifdef DEBUG
                std::printf("Initial guess for R failed (df/dr < 0 && f < 0), reducing to %E sigma\n\n", r / sigma_map[ij]);
            #endif
        }
        else{
            r = next_r;
        }
        f = get_R_rootfunc(ij, T, g, b, r);
        dfdr = get_R_rootfunc_derivative(ij, T, g, b, r);
        next_r = r - f / dfdr;
    }
    #ifdef DEBUG
        std::printf("For b = %E sigma, g = %E\n", b / sigma_map[ij], g);
        std::printf("Found R at %E sigma\n\n", next_r / sigma_map[ij]);
    #endif
    return next_r * (1 + 1e-5);
}

double KineticGas::chi(int ij, double T, double g, double b){
    const double b_cutoff = 10.0 * sigma_map[ij];
    double R_lim = get_R(ij, T, g, b_cutoff);
    double R = get_R(ij, T, g, b);
    int N_gridpoints = 50;
    double theta_lim = theta(ij, T, R_lim, g, b_cutoff, N_gridpoints);
    double theta_val = theta(ij, T, R, g, b, N_gridpoints);

    double val;
    if (abs(theta_val - theta_lim) < FLTEPS) val = 0.0;
    else val = PI - 2.0 * (theta_val - theta_lim + (PI / 2.0));
    #ifdef DEBUG
        std::printf("Limiting theta : %E pi\n", theta_lim / PI);
        std::printf("Theta : %E pi\n", theta_val / PI);
        std::printf("For b = %E sigma, g = %E\n", b / sigma_map[ij], g);
        std::printf("Computed chi = %E pi \n\n", val);
    #endif
    return val;
}
#pragma endregion

#pragma region // Bindings

#ifndef DEBUG
PYBIND11_MODULE(KineticGas_r, handle){
#else
PYBIND11_MODULE(KineticGas_d, handle){
#endif
    handle.doc() = "Is this documentation? This is documentation.";
    handle.def("cpp_tests", &cpp_tests);
    handle.def("ipow", &ipow);
    handle.def("logspace", &logspace);
    handle.def("erfspace", &erfspace);
    
    py::class_<Product>(handle, "Product")
        .def(py::init<int>())
        .def(py::init<double>())
        .def(py::init<Fac>())
        .def("eval", &Product::eval);

    py::class_<Fac>(handle, "Fac")
        .def(py::init<int>())
        .def("eval", &Fac::eval);
    
    py::class_<KineticGas>(handle, "cpp_KineticGas")
        .def(py::init<
                        std::vector<double>, 
                        std::vector<std::vector<double>>,
                        std::vector<std::vector<double>>,
                        std::vector<std::vector<double>>,
                        std::vector<std::vector<double>>,
                        int
                    >()
            )
        .def("get_A_matrix", &KineticGas::get_A_matrix)
        .def("get_delta_vector", &KineticGas::get_delta_vector)
        .def("get_reduced_A_matrix", &KineticGas::get_reduced_A_matrix)
        .def("get_alpha_vector", &KineticGas::get_alpha_vector)
        .def("A", &KineticGas::A)
        .def("A_prime", &KineticGas::A_prime)
        .def("A_trippleprime", &KineticGas::A_trippleprime)
        .def("H_ij", &KineticGas::H_ij)
        .def("H_i", &KineticGas::H_i)
        .def("H_simple", &KineticGas::H_simple)

        .def("chi", &KineticGas::chi)
        .def("get_R", &KineticGas::get_R)
        .def("potential", &KineticGas::potential)
        .def("potential_derivative_r", &KineticGas::potential_derivative_r)
        .def("potential_dblderivative_rr", &KineticGas::potential_dblderivative_rr)
        .def("omega", &KineticGas::omega)

        .def("get_R_rootfunc", &KineticGas::get_R_rootfunc)
        .def("get_R_rootfunc_derivative", &KineticGas::get_R_rootfunc_derivative)

        .def("theta", &KineticGas::theta)
        .def("theta_integrand", &KineticGas::theta_integrand)
        .def("theta_integrand_dblderivative", &KineticGas::theta_integrand_dblderivative);
}

#pragma endregion
