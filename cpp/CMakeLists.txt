cmake_minimum_required(VERSION 3.4)
set(PYBIND11_PYTHON_VERSION "3.9")
set(CMAKE_CXX_COMPILER /usr/bin/g++)
set(CMAKE_C_COMPILER /usr/bin/gcc)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
add_subdirectory(pybind11)

set(SOURCES KineticGas.cpp potentials.cpp collision_integrals.cpp Integration/Integration.cpp Factorial.cpp bindings.cpp)

if (CMAKE_BUILD_TYPE MATCHES Release)
    project(KineticGas_r)
    set(TARGET KineticGas_r)
else()
    project(KineticGas_d)
    set(TARGET KineticGas_d)
endif()

pybind11_add_module(${TARGET} ${SOURCES})