cmake_minimum_required(VERSION 3.4)
project(Cpp_KineticGas)
add_subdirectory(pybind11)
pybind11_add_module(KineticGas STATIC KineticGas.cpp Factorial.cpp)