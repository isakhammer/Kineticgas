cmake_minimum_required(VERSION 3.4)
set(PYBIND11_PYTHON_VERSION "3.9")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wno-unused-parameter -Wno-unused-const-variable")
add_subdirectory(../pybind11 release)

set(SOURCES Integration.cpp integrator_tests.cpp bindings.cpp)

if (CMAKE_BUILD_TYPE MATCHES Release)
    message("Build type is Release")
    project(Integration_r)
    set(TARGET Integration_r)
    #pybind11_add_module(Integration_r Integration.cpp integrator_tests.cpp)
else()
    message("Build type is Debug")
    project(Integration_d)
    set(TARGET Integration_d)
    #pybind11_add_module(Integration_d Integration.cpp integrator_tests.cpp)
endif()

pybind11_add_module(${TARGET} ${SOURCES})

