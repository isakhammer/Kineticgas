#pragma once
#include <vector>
#include <algorithm>
#include <memory>
#include <utility>
#include <map>
#include <iostream>

class Point{
    public:
    double x, y, z;
    Point(double x, double y) : x{x}, y{y}, z{0} {};
    Point(double x, double y, double z) : x{x}, y{y}, z{z} {};

    Point operator+(const Point& rhs){
        return Point{x + rhs.x, y + rhs.y, z + rhs.z};
    }

    void operator+=(const Point& rhs){
        x += rhs.x; y += rhs.y; z += rhs.z;
    }
};

struct Plane{
    const double A, B, C;
    Plane(double A, double B, double C) : A{A}, B{B}, C{C} {};
};

struct Line{
    const double a, b;
    Line(double a, double b) : a{a}, b{b} {};
};

Plane get_plane(const Point& p1, const Point& p2, const Point& p3); // Return vector of [A, B, C] where these are the coefficients of z = Ax + By + C
Line get_line(const Point& p1, const Point& p2); // Return vector of [A, B] where y = Ax + B
double integrate_plane_py(const Point& p1, const Point& p2, const Point& p3); // Integrate the plane interpolating (p1, p2, p3) with the triangle spanned by (p1, p2, p3)
double integrate_plane(std::shared_ptr<const Point> p1, std::shared_ptr<const Point> p2, std::shared_ptr<const Point> p3);
double eval_function(std::shared_ptr<Point> p, const int& Nx, const int& Ny,
                        double (*func)(double, double),
                        std::map<std::pair<int, int>, const double>& evaluated_points);
double eval_function(Point p, const int& Nx, const int& Ny,
                        double (*func)(double, double),
                        std::map<std::pair<int, int>, const double>& evaluated_points);
void integration_step(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2, std::shared_ptr<Point>& p3,
                        int& Nx, int& Ny, double& integral,
                        const double& dx, const double& dy,
                        int& Nxsteps, const int& Nysteps,
                        const double subdomain_dblder_limit,
                        std::map<std::pair<int, int>, const double>& evaluated_points,
                        double (*func)(double, double), std::vector<Point>& points);
double integrate_adaptive(const Point& origin,
                            const int& Nx_origin, const int& Ny_origin,
                            const int& Nx_end, const int& Ny_end,
                            const double& dx, const double& dy,
                            int& Nxsteps, const int& Nysteps,
                            const double& subdomain_dblder_limit,
                            std::map<std::pair<int, int>, const double>& evaluated_points,
                            double (*func)(double, double), std::vector<Point>& points);
double integrate2d(const Point& origin, const Point& end,
                    const double& dx, const double& dy,
                    const int& refinement_levels,
                    const double& subdomain_dblder_limit,
                    double (*func)(double, double));
std::vector<std::vector<double>> mesh2d(const Point& origin, const Point& end,
                                        const double& dx, const double& dy,
                                        const int& refinement_levels,
                                        const double& subdomain_dblder_limit,
                                        double (*func)(double, double));
double testfun(double x, double y);
double integrator_test(double origin_x, double origin_y, double end_x, double end_y,
                       double dx, double dy, int refinement_levels, double subdomain_dblder_limit);
std::vector<std::vector<double>> mesh_test(double origin_x, double origin_y, double end_x, double end_y,
                                            double dx, double dy, int refinement_levels, double subdomain_dblder_limit);
