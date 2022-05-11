# KineticGas
Implementation of Chapman-Enskog solutions to the Boltzmann equations for diffusion, thermal diffusion and conductivity for hard spheres and Mie-fluids. The file [theory.pdf](https://github.com/vegardjervell/Kineticgas/blob/main/theory.pdf) contains an excerpt of the thesis this package was created to produce, that outlines the elements of kinetic gas theory used in the package and some notable results regarding the stability, accuracy and reliability of the solutions. See the documentation for SAFT-VR-Mie at [ThermoPack](https://github.com/SINTEF/thermopack) for more details on mixing rules and the Barker-Henderson diameter.

## Dependencies
C++ module uses the [pybind11](https://github.com/pybind/pybind11) package to expose itself to the Python wrapper, removing this dependency does not amount to more than removing `bindings.cpp` from the `SOURCES` listed in `cpp/CMakeLists.txt`.

The Python wrapper requires the [ThermoPack](https://github.com/SINTEF/thermopack) python module (pyctp) and associated dependencies. The ThermoPack module is only used as a database for Mie-parameters. Removing the appropriate import statements and associated function calls will not break the code, but require that Mie-parameters are explicitly supplied.

## Setup
The package that can be installed with `pip` comes with the compiled files `KineticGas_r.so` and `KineticGas_d.so`, compiled on MacOS 10.14.6 for Python 3.9.

Build for Python 3.9 by running `bash cpp/build_kingas.sh` from the top-level directory. For a debug build, run `bash cpp/build_kingas.sh --Debug` The same script works for Linux, possibly with minor modifications. To build for different Python versions, edit the variable `PYBIND11_PYTHON_VERSION` in `cpp/CMakeLists.txt`.
For Windows, may God be with you.

Install with `pip` by running `pip install pykingas/` from the top-level directory after activating your python-installation of choice.

The Integration and meshing module can be built separately by running `bash cpp/build_integration.sh`. The integration build script also accepts the `--Debug` flag.

## Usage
Initialize a KineticGas object with the desired components and potential model, compute diffusion coefficients, thermal diffusion coefficients and thermal conductivity with the respective functions in `py_KineticGas.py`

If the argument `-d` or `-debug` or `-Debug` is passed to a script using the `pykingas` module as `python <MyScript.py> -d` the Debug build of the C++ module will be called by the python wrapper.

Some example usage can be found in the various `*_unittests.py` files.

## Acknowledgments and sources
This implementation of the Enskog solutions presented by Chapman and Cowling (*The mathematical theory of non-uniform gases* 2nd ed. Cambridge University Press, 1964) utilises the explicit summational expressions for the required bracket integrals published by Tompson, Tipton and Loyalka in *Chapman–Enskog solutions to arbitrary order in Sonine polynomials IV: Summational expressions for the diffusion- and thermal conductivity-related bracket integrals*, [European Journal of Mechanics - B/Fluids, **28**, 6, pp. 695 - 721, 2009](https://doi.org/10.1016/j.euromechflu.2009.05.002).
For a summary of the relevant theory, see the [Theory](https://github.com/vegardjervell/Kineticgas/blob/main/theory.pdf) docs.
