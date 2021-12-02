#include "KineticGas.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include <vector>
#include <algorithm>
#include <thread>

constexpr double BOLTZMANN = 1.38064852e-23;
constexpr double GAS_CONSTANT = 8.31446261815324;
constexpr double PI = 3.14159265359;
constexpr double FLTEPS = 1e-10;
#define pprintf(flt) std::printf("%f", flt); std::printf("\n")
#define pprinti(i) std::printf("%i", i); std::printf("\n")

namespace py = pybind11;

#pragma region // Helper functions

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


#pragma region // Tests
int cpp_tests(){
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
    std::vector<double> x {0.3, 0.7};
    KineticGas k{Mm, sigmaij, x, 300, 1e5, 5};
    std::vector<double> tsts {k.x1 - 0.3, k.x2 - 0.7, k.m0 - 15.6, k.sigma1 - 1.5, k.sigma2 - 2.5, k.sigma12 - 2.0};
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

#pragma region // KineticGas

#pragma region // Helper functions and initializer

int delta(int i, int j){
    if (i == j){
        return 1;
    }
    return 0;
}

double w(int l, int r){
    int f = Fac(r + 1).eval();
    if (l % 2 == 0){
        return 0.25 * (2 - ((1.0 / (l + 1)) * 2)) * f;
    }
    return 0.5 * f;
}

std::vector<std::vector<double>> KineticGas::get_A(){
    return A_matrix;
}
std::vector<double> KineticGas::get_delta(){
    return delta_vector;
}

KineticGas::KineticGas(std::vector<double> init_mole_weights, 
        std::vector<std::vector<double>> init_sigmaij, 
        std::vector<double> init_mole_fracs,
        double init_T,
        double init_p, // Pa
        int init_N) 
        : mole_weights{init_mole_weights},
        sigmaij{init_sigmaij},
        mole_fracs{init_mole_fracs},
        T{init_T},
        N{init_N},
        m0{0.0},
        A_matrix(2*init_N + 1, std::vector<double>(2 * init_N + 1)),
        delta_vector(2 * init_N + 1)
    {
    m0 = 0.0;
    for (int i = 0; i < sigmaij.size(); i++){
        sigma.push_back(sigmaij[i][i]);
        m0 += mole_weights[i];
    }
    n = init_p / (BOLTZMANN * init_T);
    x1 = mole_fracs[0];
    x2 = mole_fracs[1];
    sigma1 = sigma[0];
    sigma2 = sigma[1];
    sigma12 = sigmaij[0][1];
    M1 = mole_weights[0] / m0;
    M2 = mole_weights[1] / m0;
    
    for (int p = - N; p <= N; p++){
        for (int q = - N; q <= p; q++){
            A_matrix[p + N][q + N] = a(p, q);
            A_matrix[q + N][p + N] = A_matrix[p + N][q + N]; // Matrix is symmetric
        }
    }
    
    double delta_0 = (3.0 / (n * 2.0)) * sqrt(2 * BOLTZMANN * T / m0);
    delta_vector[N] = delta_0;
}

double KineticGas::omega(int ij, int l, int r){
    double val;
    if (ij == 1 || ij == 2){
        return pow(sigma[ij - 1], 2) * sqrt((PI * BOLTZMANN * T) / mole_weights[ij - 1]) * w(l, r);
    }
    return 0.5 * pow(sigma12, 2) * sqrt(2 * PI * BOLTZMANN * T / (m0 * M1 * M2)) * w(l, r);
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
#pragma endregion

#pragma region // Bindings
PYBIND11_MODULE(KineticGas, handle){
    handle.doc() = "This is documentation";
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
        .def(py::init<std::vector<double>, 
                        std::vector<std::vector<double>>, 
                        std::vector<double>,
                        double,
                        double,
                        int>())
        .def_property_readonly("A_matrix", &KineticGas::get_A)
        .def_property_readonly("delta_vector", &KineticGas::get_delta)
        .def("A", &KineticGas::A)
        .def("A_prime", &KineticGas::A_prime)
        .def("A_trippleprime", &KineticGas::A_trippleprime)
        .def("H_ij", &KineticGas::H_ij)
        .def("H_i", &KineticGas::H_i)
        .def("H_simple", &KineticGas::H_simple);
}

#pragma endregion