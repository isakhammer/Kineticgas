#pragma once
#include <vector>

struct Point{
    const double x, y, z;
    Point(double x, double y) : x{x}, y{y}, z{0} {};
    Point(double x, double y, double z) : x{x}, y{y}, z{z} {};
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
double integrate(const Point& p1, const Point& p2, const Point& p3); // Integrate the plane interpolating (p1, p2, p3) with the triangle spanned by (p1, p2, p3)

