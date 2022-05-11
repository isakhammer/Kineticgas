# Factorial module
Defines the types `Frac`, `Prod` and `Fac`, representing fractions, products and factorials respectively. This is required to avoid overflow and float truncation issues when evaluating fractions containing factorials that largely cancel, but where the numerator and denominator cannot be evaluated separately.

# Kineticgas module
Contains the `KineticGas` class. This class is initialized for a given binary mixture, and defines functions that solve the neccesary equations to produce the vectors and matrices required to determine diffusion coefficients and conductivities of dilute hard-sphere mixtures.

# Integration module
Functional module containing the primary functions `integrate2d()` and `mesh2d()`. See comments in code for details on usage, mathematical description still in progress but will be found in [theory.pdf](https://github.com/vegardjervell/Kineticgas/blob/main/theory.pdf).