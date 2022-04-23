#include "KineticGas.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>

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

#pragma endregion

#pragma region // Tests
int cpp_tests(){
    std::printf("We runnin!");
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
    std::invoke(potential_p, this, ij, r, theta);
}

double KineticGas::potential_derivative_r(int ij, double r, double theta){
    std::invoke(p_potential_derivative_r, this, ij, r, theta);
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
    return C_map[ij] * eps_map[ij] * ((- lr_map[ij] * pow(sigma_map[ij], lr_map[ij]) / pow(r, lr_map[ij] - 1)) + la_map[ij] * pow(sigma_map[ij], la_map[ij]) / pow(r, la_map[ij] - 1) );
}

#pragma endregion

#pragma region // Helper funcions for computing dimentionless collision integrals

double KineticGas::theta(int ij, double T, double r_prime, double g, double b){
    double lower_limit = r_prime;
    double upper_limit = lower_limit;
    double theta_0 = theta_integrand(ij, T, lower_limit, g, b);
    constexpr int N_gridpoints = 100;
    do{
        upper_limit += lower_limit;
    } while (theta_integrand(ij, T, upper_limit, g, b) > 1e-3 * theta_0);

    std::vector<double> r_grid(N_gridpoints);
    double dr = (upper_limit - lower_limit) / N_gridpoints;
    // Making a logarithmic grid
    // List is inverted when going from lin to log (so that numbers are more closely spaced at the start)
    // Therefore: Count "backwards" to invert the list again so that the smallest r is at r_grid[0]
    for (int i = 0; i < N_gridpoints; i++){
        double r = upper_limit - dr * i; // Counting backwards
        r_grid[i] = (log(r) - log(upper_limit)) * ((lower_limit - upper_limit) / (log(upper_limit) - log(lower_limit))) + lower_limit;
    }
    // Trapezoid rule (piecewise linear interpolation) integration
    double integral = 0;
    double integrand1;
    double integrand2;
    double A_coeff, B_coeff; // Coefficients for y = Ax + B, that interpolates the integrand in points (r1, r2)
    double r1 = r_grid[0];
    double r2 = r_grid[1];
    integrand1 = theta_integrand(ij, T, r1, g, b);
    integrand2 = theta_integrand(ij, T, r2, g, b);
    A_coeff = (integrand2 - integrand1) / (r2 - r1);
    B_coeff = integrand1 - A_coeff * r1;
    integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1) ; 
    for (int i = 2; i < N_gridpoints; i++){
        r1 = r_grid[i - 1];
        r2 = r_grid[i];
        integrand1 = integrand2;
        integrand2 = theta_integrand(ij, T, r2, g, b);
        A_coeff = (integrand2 - integrand1) / (r2 - r1);
        B_coeff = integrand1 - A_coeff * r1;
        integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1);
    }
    return integral;
}

double KineticGas::theta_integrand(int ij, double T, double r, double g, double b){
    // Passing dummy value "1.0" to potential. Mie potential is spherical, so not a function of last parameter (theta).
    return pow((pow(r, 4) / pow(b, 2)) * (1.0 - potential(ij, r, 1.0) / (BOLTZMANN * T * pow(g, 2))) - pow(r, 2), -0.5);
}

double KineticGas::get_R_rootfunc(int ij, double T, double g, double b, double& r){
    return (potential(ij, r, 1.0) / (BOLTZMANN * T * pow(g, 2))) + pow(b / r, 2) - 1;
}

double KineticGas::get_R_rootfunc_derivative(int ij, double T, double g, double b, double& r){
    return (potential_derivative_r(ij, r, 1.0) / (BOLTZMANN * T * pow(g, 2))) - 2 * pow(b, 2) / pow(r, 3);
}

double KineticGas::get_R(int ij, double T, double g, double b){
    // Newtons method
    double delta{100.0};
    double tol = 1e-12;
    double r = b;
    double next_r;
    do{
        next_r = r - get_R_rootfunc(ij, T, g, b, r) / get_R_rootfunc_derivative(ij, T, g, b, r);
        delta = abs(r - next_r);
        r = next_r;
    }
    while (delta > tol);
    return r * (1 + 1e-5);
}

double KineticGas::chi(int ij, double T, double g, double b){
    double R = get_R(ij, T, g, b);
    return PI - 2 * theta(ij, T, R, g, b);
}
#pragma endregion

#pragma region // Bindings
PYBIND11_MODULE(KineticGas, handle){
    handle.doc() = "Is this documentation? This is documentation.";
    handle.def("cpp_tests", &cpp_tests);
    handle.def("ipow", &ipow);
    
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
        
        .def("theta", &KineticGas::theta)
        .def("chi", &KineticGas::chi)
        .def("get_R", &KineticGas::get_R);
}

#pragma endregion
