"""
    PeriodicInterpolation1D

Module for 1D periodic interpolation on uniform grids using various methods.

This module provides different interpolation techniques for periodic 1D domains:
- **Lagrange**: Global FFT-based Lagrange polynomials (spectral accuracy)
- **BSpline**: B-spline basis functions (smooth, controllable accuracy)
- **Spectral**: Pure Fourier spectral interpolation (exponential accuracy)
- **FastLagrange**: Local Lagrange polynomials (fast, limited range)

All methods assume periodic boundary conditions on the input data.

## Authors
- Original Fortran: Klaus Reuter (MPCDF), Katharina Kormann (RUB), Michel Mehrenberger (I2M)
- Julia Translation: Pierre Navaro (IRMAR)

## Reference
Based on formulas from Abramowitz and Stegun: Handbook of Mathematical Functions, Chapter 25.2

## Quick Start

```julia
using PeriodicInterpolation1D

# Create sample data on 100-point periodic grid
n = 100
x = 2π .* (0:n-1) ./ n
f = sin.(x)
f_interp = zeros(n)

# Interpolate using Lagrange polynomials
lagr = Lagrange(n, 5)
interpolate!(f_interp, lagr, f, 0.5)  # shift by 0.5 grid points

# Or use B-splines
bspl = BSpline(n, 4)
interpolate!(f_interp, bspl, f, 0.5)

# Or use spectral method
spec = Spectral(n)
interpolate!(f_interp, spec, f, 0.5)
```

## Method Comparison

| Method | Accuracy | Speed | Best Use |
|--------|----------|-------|----------|
| Lagrange | Spectral | Medium | Smooth functions, high accuracy needed |
| BSpline | Polynomial | Medium | Smooth interpolation, controllable accuracy |
| Spectral | Exponential | Medium | Very smooth data, maximum accuracy |
| FastLagrange | Polynomial | Fast | Small shifts, performance critical |
"""
module PeriodicInterpolation1D

using DocStringExtensions
using FFTW

# Export main interpolation objects and function
export Lagrange, BSpline, Spectral, FastLagrange
export interpolate!

include("bsplines.jl")
include("spectral.jl")
include("fast_lagrange.jl")
include("lagrange.jl")
include("cubic_interp.jl")

end # module
