/*
Contains functions for computing the collision integrals
Common variables are:
    ij : collision type (1 for like molecules of type 1, 2 for type 2, and 12 or 21 for collisions of unlike molecules
    T : Temperature [K]
    l, r : Collision integral indices (in omega and w_<potential_type> functions)
    g : Dimentionless relative velocity of molecules
    b : "Impact paramter", closest distance of approach of the centers of mass if molecules were non-interacting point particles
    theta : Angular polar coordinate. Theta = 0 when the particles are infinitely far away from each other, before the collision,
        theta(t1) is the angle between the line between the particles before interaction begins (at infinite distance), and the line between the particles at t = t1
    r : Radial polar coordinate. Distance from the center of mass of one particle to the other.
    R : Actual distance of closest approach. Corresponds to dr/dtheta = 0
    chi : Deflection angle, corresponds to pi - 2 * theta(R)
*/

#include "KineticGas.h"
#include "Integration/Integration.h"

#pragma region // Collision integrals for various potentials

double KineticGas::omega(int ij, int l, int r){
    std::printf("Computing omega for (%i, %i)\n", l, r);
    double w = std::invoke(w_p, this, ij, T, l, r); // w_p is a pointer to the dimentionless collision integral corresponding to this.potential_mode
    if (ij == 1 || ij == 2){
        return pow(sigma[ij - 1], 2) * sqrt((PI * BOLTZMANN * T) / mole_weights[ij - 1]) * w;
    }
    return 0.5 * pow(sigma12, 2) * sqrt(2 * PI * BOLTZMANN * T / (m0 * M1 * M2)) * w;
}

// Dimentionless collision integral for a Hard-sphere potential (analytic)
double KineticGas::w_HS(int ij, double T, int l, int r){
    int f = Fac(r + 1).eval();
    if (l % 2 == 0){
        return 0.25 * (2 - ((1.0 / (l + 1)) * 2)) * f;
    }
    return 0.5 * f;
}

// Dimentionless collision integral for a Mie-potential
double KineticGas::w_spherical(int ij, double T, int l, int r){
    Point origin{1e-5, 1e-5};
    Point end{5, 3};
    double dx{0.05}, dy{0.05};
    int refinement_levels{4};
    double subdomain_dblder_limit{0.05};
    std::function<double(int, double, double, double, int, int)> f = std::bind(&KineticGas::w_spherical_integrand, this,
                                                                                std::placeholders::_1, std::placeholders::_2, 
                                                                                std::placeholders::_3, std::placeholders::_4,
                                                                                std::placeholders::_5, std::placeholders::_6);
    return integrate2d(origin, end,
                        dx, dy,
                        refinement_levels,
                        subdomain_dblder_limit,
                        ij, T, l, r,
                        f);
}

double KineticGas::w_spherical_integrand(const int& ij, const double& T, 
                                        const double& g, const double& b,
                                        const int& l, const int& r){ // Using b = b / sigma to better scale the axes. Multiply the final integral by sigma.
    const double chi_val = chi(ij, T, g, b * sigma_map[ij]);
    return 2 * exp(- pow(g, 2)) * pow(g, 2.0 * r + 3.0) * (1 - pow(cos(chi_val), l)) * b;
};

#pragma region // Helper funcions for computing dimentionless collision integrals

double KineticGas::theta(const int ij, const double T, const double g, const double b){
    // Compute deflection angle for a collision
    if (b / sigma_map[ij] > 10) return PI / 2;
    if (b / sigma_map[ij] < 1e-3) return 0;
    double R = get_R(ij, T, g, b);
    return theta_integral(ij, T, R, g, b) - theta_lim(ij, T, g) + PI / 2;
}

double KineticGas::theta_lim(const int ij, const double T, const double g){
    double b = 10 * sigma_map[ij];
    double R = get_R(ij, T, g, b);
    return theta_integral(ij, T, R, g, b);
}

double KineticGas::theta_integral(const int ij, const double T, const double R, const double g, const double b){

    double lower_cutoff_factor = 1e-7;
    double r_prime = (1 + lower_cutoff_factor) * R;

    double step_eps_tol = 1e-6; // Absolute Error tolerance in each integration step
    double rel_eps_tol = 1e-3; // Relative error tolerance for total integration
    double abs_eps_tot; // Total absolute error

    double total_integral;
    double integral_09; // Integral at 90% completion
    double r09;
    double convergence_threshold = 0.99; // Integral has converged when integral_09 >= convergence_threshold * total_integral
    double A_coeff, B_coeff; // Coefficients for y = Ax + B, that interpolates the integrand in points (r1, r2)

    int N_total_integration_steps;
    double r2;
    double r1; // will be overwritten by r2 at the start of the first iteration
    double delta_r;
    double integrand2;
    double integrand1;

    // Trapezoid rule (piecewise linear interpolation) integration
    do { // Compute integral, if total relative error exceeds tolerance, reduce stepwise tolerance and recompute.
        total_integral = 0;
        abs_eps_tot = 0;
        N_total_integration_steps = 0;
        r2 = r_prime;
        integrand2 = theta_integrand(ij, T, r2, g, b);
        for (; N_total_integration_steps < 100; N_total_integration_steps++){ // Start by doing 100 integration steps
            r1 = r2;
            double d3tdr3 = theta_integrand_dblderivative(ij, T, r1, g, b);
            double M = max(abs(d3tdr3), pow(0.1 * sigma_map[ij], 3) / (12.0 * step_eps_tol)); // Maximum step length is 0.1 * sigma
            delta_r = pow(12.0 * step_eps_tol / M, 1.0 / 3.0); // Ensure error is less than tolerance in each step
            r2 = r1 + delta_r;

            integrand1 = integrand2;
            integrand2 = theta_integrand(ij, T, r2, g, b);
            A_coeff = (integrand2 - integrand1) / (r2 - r1);
            B_coeff = integrand1 - A_coeff * r1;
            total_integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1);
            abs_eps_tot += theta_integrand_dblderivative(ij, T, r1, g, b) * pow(r2 - r1, 3) / 12.0;

        }
        int N_increment_points;
        do{ // Compute integral, continue until convergence
            r09 = r2;
            integral_09 = total_integral;
            N_increment_points = (int) (N_total_integration_steps / 9.0) + 0.5;

            for (int i = 0; i < N_increment_points; i++){
                r1 = r2;
                double d3tdr3 = theta_integrand_dblderivative(ij, T, r1, g, b);
                double M = max(abs(d3tdr3), pow(0.1 * sigma_map[ij], 3) / (12.0 * step_eps_tol)); // Maximum step length is 0.1 * sigma
                delta_r = pow(12.0 * step_eps_tol / M, 1.0 / 3.0); // Ensure error is less than tolerance in each step
                r2 = r1 + delta_r;

                integrand1 = integrand2;
                integrand2 = theta_integrand(ij, T, r2, g, b);
                A_coeff = (integrand2 - integrand1) / (r2 - r1);
                B_coeff = integrand1 - A_coeff * r1;
                total_integral += A_coeff * (pow(r2, 2) - pow(r1, 2)) / 2 + B_coeff * (r2 - r1);
                abs_eps_tot += theta_integrand_dblderivative(ij, T, r1, g, b) * pow(r2 - r1, 3) / 12.0;

                N_total_integration_steps++;
            }

            #ifdef DEBUG
                if ((integral_09 < convergence_threshold * total_integral) || (isnan(total_integral))){
                    std::printf("theta = %E pi\n", total_integral / PI);
                    std::printf("r = %E sigma\n", r2 / sigma_map[ij]);
                    std::printf("r - r09 = %E sigma\n", (r2 - r09) / sigma_map[ij]);
                    std::printf("Change in final 10%% is : %E %% of final value\n", (1 - integral_09 / total_integral) * 100.0);
                    std::printf("N_gridpoints, increment, and increment %% are %i, %i, %E %%\n", N_total_integration_steps, N_increment_points, (100.0 * N_increment_points) / N_total_integration_steps);
                    std::printf("Relative error is : %E %%\n\n", 100.0 * abs_eps_tot / total_integral);
                }
            #endif

        } while (integral_09 < convergence_threshold * total_integral);

        if ((isnan(total_integral) || isinf(total_integral))){
            if (lower_cutoff_factor > 1e-2) throw "Something wrong...";
            lower_cutoff_factor *= 10;
            r_prime = (1 + lower_cutoff_factor) * R;
            total_integral = 1;
            abs_eps_tot = rel_eps_tol + 1; // Ensure that the loop continues
        }
        else{
            step_eps_tol *= 0.5;
        }
        #ifdef DEBUG
            if (abs_eps_tot / total_integral > rel_eps_tol){
                std::printf("Relative error is less than %E, Tolerance is %E\n", abs_eps_tot / total_integral, rel_eps_tol);
                std::printf("Reducing step tolerance to %E\n", step_eps_tol);
                std::printf("Lower cutoff factor is %E\n\n", lower_cutoff_factor);
            }
        #endif

    } while (abs_eps_tot / total_integral > rel_eps_tol || step_eps_tol < 1e-10);

    if (step_eps_tol < 1e-10){
        std::printf("\nWarning : theta integral could not achieve expected precicion\n");
        std::printf("Upper limit for relative error is %E, tolerance is %E\n", abs_eps_tot / total_integral, rel_eps_tol);
        std::printf("Integral parameters are :\nij = %i \nT = %E \nb = %E sigma\ng = %E \nR = %E sigma\n\n", ij, T, b / sigma_map[ij], g, r_prime / sigma_map[ij]);
    }

    #ifdef DEBUG
        std::printf("For b = %E sigma, g = %E\n", b / sigma_map[ij], g);
        std::printf("init, start, end (sigma) = %E, %E, %E\n", R / sigma_map[ij], r_prime / sigma_map[ij], r2 / sigma_map[ij]);
        std::printf("Total gridpoints = %i\n", N_total_integration_steps);
        std::printf("Computed theta = %E pi\n", total_integral / PI);
        std::printf("Relative error is less than %E %%, tolerance is %E %% \n\n", 100 * abs_eps_tot / total_integral, 100 * rel_eps_tol);
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

    double val = (3.0 / 4.0) * pow(core, -2.5) * pow(core_prime, 2) - 0.5 * pow(core, - 1.5) * core_dblprime;
    #ifdef DEBUG
        if (val < 0){
            std::printf("\nd3tdr3 at r = %E sigma\n", r / sigma_map[ij]);
            std::printf("val = %E\n\n", val);
        }
    #endif
    return val;
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
    return next_r;
}

double KineticGas::chi(int ij, double T, double g, double b){
    if (b / sigma_map[ij] > 10) return 0;
    double t = theta(ij, T, g, b);
    double val = PI - 2.0 * t;
    #ifdef DEBUG
        std::printf("For b = %E sigma, g = %E\n", b / sigma_map[ij], g);
        std::printf("Computed chi = %E pi \n\n", val);
    #endif
    return val;
}

double KineticGas::chi_HS(int ij, double T, double g, double b){
    if (b >= sigma_map[ij]) return 0;
    return acos(1 - 2 * (1 - pow(b / sigma_map[ij], 2)));
}
