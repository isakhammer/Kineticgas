/*
File containing various potential models implemented in KineticGas
Spherical potentials must supply a function, as well as the first and second derivative with respect to r
*/
#include "KineticGas.h"

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
    // To get this, start with a potential that has a second derivative f''(r) = (sigma / r)^22 + A
    // Then integrate the function and require that f''(sigma) = f'(sigma) = f(sigma) = 0
    if (r > sigma_map[ij]){
        return 0.0;
    }
    return pow(sigma_map[ij] / r, 20) - (20.0 * 21.0 / 2) * pow(r / sigma_map[ij], 2) + 20.0 * 22.0 * (r / sigma_map[ij]) + 20.0 * ((21.0 / 2.0) - 22.0) - 1.0; // Force continiuous function
}

double KineticGas::HS_potential_derivative(int ij, double r, double theta){
    if (r > sigma_map[ij]){
        return 0.0;
    }
    return - 20.0 * pow(sigma_map[ij], 20) / pow(r, 21) - 20.0 * 21.0 * r / pow(sigma_map[ij], 2) + 20.0 * 22.0 / sigma_map[ij]; // Force continiuous first derivative
}

double KineticGas::HS_potential_dblderivative_rr(int ij, double r, double theta){
    if (r > sigma_map[ij]){
        return 0.0;
    }
    return 20.0 * 21.0 * pow(sigma_map[ij], 20) / pow(r, 22) - 20.0 * 21.0 / pow(sigma_map[ij], 2); // Force continiuous second derivative
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