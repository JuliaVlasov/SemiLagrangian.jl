# Getting Started

Quick start guide for using `PeriodicInterpolation1D`.

## Installation

Add the package to your Julia environment:

```julia
julia> using Pkg
julia> Pkg.add(url=\"https://github.com/JuliaVlasov/PeriodicInterpolation1D.jl\")
```

Or in the REPL in `pkg` mode:

```
pkg> add https://github.com/JuliaVlasov/PeriodicInterpolation1D.jl
```

## Basic Usage

### 1. Simple Example with Lagrange Interpolation

```julia
using PeriodicInterpolation1D
using Plots  # optional, for visualization

# Create a periodic grid
n = 100
x = 2π .* (0:n-1) ./ n

# Create a smooth function
f = sin.(x)

# Create output array
f_interp = zeros(n)

# Create interpolation object with 5-point stencil (order-4 polynomial)
lagr = Lagrange(n, 5)

# Interpolate with shift of 0.5 grid points
interpolate!(f_interp, lagr, f, 0.5)

println(\"Interpolated values computed: \", length(f_interp), \" points\")
```

### 2. Comparing Different Methods

Interpolate the same function using all four methods:

```julia
using PeriodicInterpolation1D

n = 100
x = 2π .* (0:n-1) ./ n
f = sin.(x)
shift = 0.3

# Output vectors
f_lagr = zeros(n)
f_bspl = zeros(n)
f_spec = zeros(n)
f_fast = zeros(n)

# Create interpolators
lagr = Lagrange(n, 7)
bspl = BSpline(n, 4)
spec = Spectral(n)
fast = FastLagrange(7)

# Perform interpolation
interpolate!(f_lagr, lagr, f, shift)
interpolate!(f_bspl, bspl, f, shift)
interpolate!(f_spec, spec, f, shift)
interpolate!(f_fast, fast, f, shift)

# Compare errors with analytical shifted solution
x_shifted = mod.(x .- shift, 2π)
f_exact = sin.(x_shifted)

err_lagr = maximum(abs.(f_lagr .- f_exact))
err_bspl = maximum(abs.(f_bspl .- f_exact))
err_spec = maximum(abs.(f_spec .- f_exact))
err_fast = maximum(abs.(f_fast .- f_exact))

println(\"Maximum absolute errors:\")
println(\"  Lagrange:    \", err_lagr)
println(\"  BSpline:     \", err_bspl)
println(\"  Spectral:    \", err_spec)
println(\"  FastLagrange:\", err_fast)
```

### 3. Multiple Interpolations

Efficiently perform many interpolations by reusing objects:

```julia
using PeriodicInterpolation1D

n = 1000
f = randn(n)  # random data

# Create interpolator once
lagr = Lagrange(n, 5)

# Perform many interpolations efficiently
shifts = range(0, 2π, length=100)
results = []

for shift in shifts
    f_interp = zeros(n)
    interpolate!(f_interp, lagr, f, shift)
    push!(results, f_interp)
end

println(\"Computed \", length(results), \" interpolations\")
```

### 4. Cosine Interpolation Example

```julia
using PeriodicInterpolation1D

# High-frequency cosine function
n = 128
x = 2π .* (0:n-1) ./ n
k = 10  # wave number
f = cos.(k .* x)

# Try different methods
f_interp = zeros(n)

# Spectral method
spec = Spectral(n)
interpolate!(f_interp, spec, f, 0.25)

# Compute exact solution
x_shifted = mod.(x .- 0.25, 2π)
f_exact = cos.(k .* x_shifted)

error = maximum(abs.(f_interp .- f_exact))
println(\"Spectral method error: \", error)
```

### 5. Large-scale Problem

```julia
using PeriodicInterpolation1D, BenchmarkTools

# Large problem
n = 10_000
f = randn(n)

# Time different methods
lagr = Lagrange(n, 5)
bspl = BSpline(n, 4)
spec = Spectral(n)
fast = FastLagrange(7)

f_out = zeros(n)
shift = 0.25

println(\"Performance comparison for n=$n, shift=$shift:\")
@time interpolate!(f_out, lagr, f, shift)
@time interpolate!(f_out, bspl, f, shift)
@time interpolate!(f_out, spec, f, shift)
@time interpolate!(f_out, fast, f, shift)
```

## Common Pitfalls

### 1. Mutable Output Array

The output array is modified **in-place**. Make sure to provide a pre-allocated array:

```julia
lagr = Lagrange(100, 5)
u = randn(100)

# ✓ Correct
u_out = zeros(100)
interpolate!(u_out, lagr, u, 0.5)

# ✗ Wrong - will not work as expected
result = interpolate!(nothing, lagr, u, 0.5)  # Error!
```

### 2. Array Size Mismatch

All arrays must have the same length:

```julia
lagr = Lagrange(100, 5)
u = randn(100)

# ✓ Correct
u_out = zeros(100)
interpolate!(u_out, lagr, u, 0.5)

# ✗ Wrong - length mismatch
u_out_wrong = zeros(50)
interpolate!(u_out_wrong, lagr, u, 0.5)  # Error!
```

### 3. FastLagrange Stability

Remember that `FastLagrange` is unstable for large shifts:

```julia
fast = FastLagrange(7)
u = randn(100)
u_out = zeros(100)

# ✓ Good - small shift
interpolate!(u_out, fast, u, 0.1)

# ⚠️ Warning - large shift, inaccurate results
interpolate!(u_out, fast, u, 10.5)  # Large shift, results unreliable!
```

### 4. BSpline Order Must Be Even

```julia
# ✓ Correct
bspl = BSpline(100, 4)  # order 4 is even

# ✗ Wrong - odd order not allowed
# bspl = BSpline(100, 5)  # error!
```

## Next Steps

- See [Interpolation Methods](@ref) for detailed method descriptions
- Check the [Reference](@ref reference) for API documentation
- Look at examples in the `examples/` directory
- Read the original Fortran documentation from [SeLaLib](https://selalib.github.io)

## Performance Considerations

1. **FFT-based vs Local Methods**
   - Spectral, Lagrange, BSpline: Use FFT, good for multiple interpolations
   - FastLagrange: Direct computation, good for single interpolations

2. **Memory Usage**
   - Spectral: ~2n complex numbers (FFT workspace)
   - Lagrange: ~2n complex numbers (FFT workspace)
   - BSpline: ~2n complex numbers (FFT workspace)
   - FastLagrange: ~stencil_size floats (minimal overhead)

3. **Accuracy vs Speed Trade-off**
   - Higher polynomial order → better accuracy but slower
   - For most problems, order 5-7 is a good balance
   - Use FastLagrange only for performance-critical small-shift cases
