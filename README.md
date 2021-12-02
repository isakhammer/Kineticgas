# Kineticgas
Implementation of Enskog solutions for diffusion, thermal diffusion and conductivity

## Dependencies
C++ module can be compiled as is.
The Python extension requires the [Thermopack](https://github.com/SINTEF/thermopack) python module (pyctp) and associated dependencies.

## Acknowledgments and sources
This implementation of the Enskog solutions presented by Chapman and Cowling (*The mathematical theory of non-uniform gases* 2nd ed. Cambridge University Press, 1964) utilises the explicit summational expressions for the required bracket integrals published by Tompson, Tipton and Loyalka in *Chapman–Enskog solutions to arbitrary order in Sonine polynomials IV: Summational expressions for the diffusion- and thermal conductivity-related bracket integrals*, [European Journal of Mechanics - B/Fluids, **28**, 6, pp. 695 - 721, 2009](https://doi.org/10.1016/j.euromechflu.2009.05.002).
