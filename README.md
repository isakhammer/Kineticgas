# Kineticgas
Implementation of Enskog solutions for diffusion, thermal diffusion and conductivity

## Dependencies
C++ module uses the [pybind11](https://github.com/pybind/pybind11) package to expose itself to the Python wrapper, removing this dependency does not ammount to more than deleting the final lines in KineticGas.cpp and removing the appropriate `include` statements.

The Python extension requires the [ThermoPack](https://github.com/SINTEF/thermopack) python module (pyctp) and associated dependencies. The ThermoPack module is only used as a database for Mie-parameters. Removing the appropriate import statements and associated function calls will not break the code, but require that Mie-parameters are explicitly supplied.

## Acknowledgments and sources
This implementation of the Enskog solutions presented by Chapman and Cowling (*The mathematical theory of non-uniform gases* 2nd ed. Cambridge University Press, 1964) utilises the explicit summational expressions for the required bracket integrals published by Tompson, Tipton and Loyalka in *Chapman–Enskog solutions to arbitrary order in Sonine polynomials IV: Summational expressions for the diffusion- and thermal conductivity-related bracket integrals*, [European Journal of Mechanics - B/Fluids, **28**, 6, pp. 695 - 721, 2009](https://doi.org/10.1016/j.euromechflu.2009.05.002).
