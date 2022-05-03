#include "KineticGas.h"
#include "Factorial.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"

namespace py = pybind11;

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