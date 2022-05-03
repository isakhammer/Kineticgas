#pragma once
#include "Integration.h"
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <map>
#include <math.h>
#include <functional>

// mesh_step and mesh_adaptive are together the exact same algorithm as integration_step and integrate_adaptive
// Defined in Integrator.h. Rather than compute the integral, mesh_adaptive returns a vector of all the integration
// points, in the order they are evaluated.
// mesh2d is the corresponding "mirror function" for integrate2d.
void mesh_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        std::function<double(double, double)> func, std::vector<Point>& points);

void mesh_adaptive(const Point& origin,
                    const int& Nx_origin, const int& Ny_origin,
                    const int& Nx_end, const int& Ny_end,
                    const double& dx, const double& dy,
                    int& Nxsteps, const int& Nysteps,
                    const double& subdomain_dblder_limit,
                    std::map<std::pair<int, int>, const double>& evaluated_points,
                    std::function<double(double, double)> func, std::vector<Point>& points);

std::vector<std::vector<double>> mesh2d(const Point& origin, const Point& end,
                                        const double& dx, const double& dy,
                                        const int& refinement_levels,
                                        const double& subdomain_dblder_limit,
                                         std::function<double(double, double)> func);

double testfun(double x, double y);
double testfun_linear(double x, double y);
double integrator_test(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit);
double integrator_test_linear(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit);
std::vector<std::vector<double>> mesh_test(double origin_x, double origin_y, double end_x, double end_y,
                                            double dx, double dy, int refinement_levels, double subdomain_dblder_limit);