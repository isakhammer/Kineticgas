#include "Integration.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include <vector>
#include <algorithm>
#include <memory>

#define FLTEPS 1e-12

namespace py = pybind11;

Plane get_plane(const Point& p1, const Point& p2, const Point& p3){
    // a, b, c are components of vector normal to plane
    const double a = - p1.y * p3.z - p2.y * p1.z + p2.y * p3.z + p1.y * p2.z + p3.y * p1.z - p3.y * p2.z;
    const double b = p1.x * p3.z + p2.x * p1.z - p2.x * p3.z - p1.x * p2.z - p3.x * p1.z + p3.x * p2.z;
    const double c = - p1.x * p3.y - p2.x * p1.y + p2.x * p3.y + p1.x * p2.y + p3.x * p1.y - p3.x * p2.y;
    const double d = - (a * p1.x + b * p1.y + c * p1.z);
    return Plane{- a / c, - b / c, - d / c};
}

Line get_line(const Point& p1, const Point& p2){
    const double A = (p2.y - p1.y) / (p2.x - p1.x);
    const double B = p2.y - A * p2.x;
    return Line{A, B};
}

double integrate(const Point& p1_in, const Point& p2_in, const Point& p3_in){
    const Point* p1 = &p1_in;
    const Point* p2 = &p2_in;
    const Point* p3 = &p3_in;
    const Plane plane = get_plane(*p1, *p2, *p3);
    if ((p1->x <= p2->x) && (p1->x <= p3->x)){ // Hvis p1_in.x er minst ...
        if (p2->x > p3->x){
            std::swap(p2, p3);
        }
    }
    else if ((p1->x >= p2->x) && (p1->x >= p3->x)){ // Hvis p1.x er størst ...
        std::swap(p1, p3); // Nå er p3 størst
        if (p1->x > p2->x){
            std::swap(p1, p2);
        }
    }
    else{ // p1_in.x er midterst
        if (p2->x < p3->x){
            std::swap(p1, p2);
        }
        else{
            std::swap(p1, p3);
            std::swap(p2, p3);
        }
    }


    std::shared_ptr<Line> l_upper_1{new Line{get_line(*p1, *p3)} }; // y    p1 * * * * * * p3   // Assuming this configuration
    std::shared_ptr<Line> l_lower_1{new Line{get_line(*p1, *p2)} }; // ^        *       *       // Swapping if y2 > y3
    std::shared_ptr<Line> l_upper_2{l_upper_1};                     // |          *   *         //
    std::shared_ptr<Line> l_lower_2{new Line{get_line(*p2, *p3)} }; //  => x     p2 *           //
    
    if (p3->y < p2->y){
        std::swap(l_lower_1, l_upper_1);
        std::swap(l_lower_2, l_upper_2);
    }
    
    double integral = 0;
    if (abs(p1->x - p2->x) > FLTEPS){ // Integral from x1 to x2

        const double A31_star = l_upper_1->a - l_lower_1->a;
        const double B31_star = l_upper_1->b - l_lower_1->b;
        const double A31_dblstar = pow(l_upper_1->a, 2) - pow(l_lower_1->a, 2);
        const double B31_dblstar = 2 * (l_upper_1->a * l_upper_1->b - l_lower_1->a * l_lower_1->b);
        const double C31_dblstar = pow(l_upper_1->b, 2) - pow(l_lower_1->b, 2);

        const double A12_tilde = plane.A * A31_star + plane.B * A31_dblstar / 2;
        const double B12_tilde = plane.A * B31_star + plane.C * A31_star + plane.B * B31_dblstar / 2;
        const double C12_tilde = plane.C * B31_star + plane.B * C31_dblstar / 2;

        integral += (A12_tilde / 3) * (pow(p2->x, 3) - pow(p1->x, 3)) + (B12_tilde / 2) * (pow(p2->x, 2) - pow(p1->x, 2)) + C12_tilde * (p2->x - p1->x);
    }
    if (abs(p2->x - p3->x) > FLTEPS){ // Integral from x2 to x3

        const double A23_star = l_upper_2->a - l_lower_2->a;
        const double B23_star = l_upper_2->b - l_lower_2->b;
        const double A23_dblstar = pow(l_upper_2->a, 2) - pow(l_lower_2->a, 2);
        const double B23_dblstar = 2 * (l_upper_2->a * l_upper_2->b - l_lower_2->a * l_lower_2->b);
        const double C23_dblstar = pow(l_upper_2->b, 2) - pow(l_lower_2->b, 2);

        const double A23_tilde = plane.A * A23_star + plane.B * A23_dblstar / 2;
        const double B23_tilde = plane.A * B23_star + plane.C * A23_star + plane.B * B23_dblstar / 2;
        const double C23_tilde = plane.C * B23_star + plane.B * C23_dblstar / 2;

        integral += (A23_tilde / 3) * (pow(p3->x, 3) - pow(p2->x, 3)) + (B23_tilde / 2) * (pow(p3->x, 2) - pow(p2->x, 2)) + C23_tilde * (p3->x - p2->x);
    }
    return integral;


}

#ifndef DEBUG
PYBIND11_MODULE(Integration_r, handle){
#else
PYBIND11_MODULE(Integration_d, handle){
#endif
    handle.doc() = "Integration module";
    handle.def("get_plane", &get_plane);
    handle.def("get_line", &get_line);
    handle.def("integrate", &integrate);

    py::class_<Point>(handle, "Point")
        .def(py::init<double, double, double>())
        .def(py::init<double, double>())
        .def_readonly("x", &Point::x)
        .def_readonly("y", &Point::y)
        .def_readonly("z", &Point::z);

    py::class_<Plane>(handle, "Plane")
        .def(py::init<double, double, double>())
        .def_readonly("A", &Plane::A)
        .def_readonly("B", &Plane::B)
        .def_readonly("C", &Plane::C);

    py::class_<Line>(handle, "Line")
        .def(py::init<double, double>())
        .def_readonly("a", &Line::a)
        .def_readonly("b", &Line::b);
}