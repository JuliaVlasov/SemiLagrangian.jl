```@meta
CurrentModule = SemiLagrangian
```

# 1D Periodic Interpolation on Uniform Grids

Periodic interpolation methods available in `SemiLagrangian`.

Author: Translated from Fortran to Julia by Pierre Navaro with [claude.ai](https://claude.ai/).
Original Authors: Klaus Reuter (MPCDF), Katharina Kormann (RUB) and Michel Mehrenberger (I2M)

## Description

Module for different 1D interpolation on a uniform grid.
This is an implementation for equidistant grids that exploits the simplifications in this
special case in order to be faster.
