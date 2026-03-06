# Interpolation Methods

This page describes the different interpolation methods available in `PeriodicInterpolation1D`.

## Overview

The module provides five interpolation methods for periodic 1D grids on uniform meshes:

1. **[Lagrange](@ref)** - Global FFT-based Lagrange polynomials
2. **[BSpline](@ref)** - B-spline basis functions  
3. **[Spectral](@ref)** - Pure Fourier spectral interpolation
4. **[FastLagrange](@ref)** - Local Lagrange polynomials
5. **[CubicSpline](@ref)** - Cubic-spline interpolation (uses `FastInterpolations.jl` for speed)

## Method Comparison

| Method | Accuracy | Speed | Memory | Best For |
|--------|----------|-------|--------|----------|
| Lagrange | Spectral | Medium | Medium | Smooth functions, high accuracy |
| BSpline | Polynomial order | Medium | Medium | Smooth interpolation, flexibility |
| Spectral | Exponential | Medium | Medium | Very smooth data, maximum accuracy |
| FastLagrange | Polynomial order | Fast | Low | Small shifts, performance-critical code |
| CubicSpline | Cubic | Fast | Medium | Smooth cubic-spline interpolation (fast with FastInterpolations.jl) |

## Detailed Method Descriptions

### Lagrange Interpolation

The `Lagrange` method uses FFT-based computation with global Lagrange polynomials:

```julia
using PeriodicInterpolation1D

n = 100
x = 2π .* (0:n-1) ./ n
f = sin.(x)
f_interp = zeros(n)

lagr = Lagrange(n, 5)  # 5-point stencil → order-4 polynomial
interpolate!(f_interp, lagr, f, 0.5)  # shift by 0.5 grid points
```

**Advantages:**
- Spectral accuracy for smooth periodic functions
- Arbitrary shift values supported
- FFT leverage for efficiency

**Disadvantages:**
- Requires FFT of entire domain
- Memory overhead for FFT arrays
- Slower for small problems

### B-Spline Interpolation

The `BSpline` method uses B-spline basis functions with even polynomial orders:

```julia
bspl = BSpline(n, 4)  # even order: cubic B-splines (C² continuity)
interpolate!(f_interp, bspl, f, 0.3)
```

**Supported Orders:**

- 2: linear B-splines (C⁰ continuity)
- 4: cubic B-splines (C² continuity)  
- 6: quintic B-splines (C⁴ continuity)
- 8: septic B-splines (C⁶ continuity)

**Advantages:**

- Smooth interpolation with controlled continuity
- Good balance between accuracy and speed
- Efficient FFT implementation

**Disadvantages:**

- Orders must be even
- Requires even order → 2nd order and higher

### Spectral Interpolation

The `Spectral` method uses pure Fourier space shifting for maximum accuracy:

```julia
spec = Spectral(n)
interpolate!(f_interp, spec, f, 0.25)
```

**Advantages:**
- Exponential accuracy for smooth periodic functions
- Simplest implementation
- Exact for band-limited functions

**Disadvantages:**
- FFT overhead
- Prone to Gibbs phenomena for non-smooth data
- Not suitable for discontinuous functions

### FastLagrange Interpolation

The `FastLagrange` method uses local stencil-based Lagrange polynomials optimized for speed:

```julia
fast_lag = FastLagrange(7)  # 7-point stencil
interpolate!(f_interp, fast_lag, f, 0.1)  # small shift recommended
```

**Supported Stencil Widths:**
- 3: order-2 polynomial (1-point half-width)
- 5: order-4 polynomial (2-point half-width)
- 7: order-6 polynomial (3-point half-width)
- 9: order-8 polynomial (4-point half-width)
- 11: order-10 polynomial (5-point half-width)

**⚠️ Stability Warning:**
- **Limited to small shifts**: Best accuracy for `|shift| ≤ 1.0` grid spacing
- **Runge's phenomenon**: Large shifts cause weight oscillation and numerical instability
- **Local approximation**: Only nearby points contribute to interpolation

**Advantages:**
- Fastest method for small shifts
- Minimal memory overhead
- No FFT required

**Disadvantages:**
- Unstable for large shifts
- Lower accuracy than spectral methods
- Must keep shifts small

### Cubic Spline Interpolation

The `CubicSpline` method performs cubic-spline interpolation on a periodic
uniform grid. Internally this adapter delegates to `FastInterpolations.jl`'s
`cubic_interp!` routine and applies periodic boundary conditions. It offers a
good balance between smoothness and performance for many practical problems.

```julia
cs = CubicSpline(n)
interpolate!(u_out, cs, u, 0.25)
```

**Advantages:**
- Smooth cubic interpolation with C^2 continuity
- Fast implementations available via `FastInterpolations.jl`

**Disadvantages:**
- Requires the `FastInterpolations.jl` dependency for the high-performance
   backend used here
- Not spectrally accurate like FFT-based methods

## Accuracy Comparison

For a smooth periodic function (e.g., `sin(2πx)`), typical accuracy (L∞ error) for various methods with 100 grid points and 0.5 grid-point shift:

```
Method         | Error magnitude
Spectral       | ~10^-13 (machine precision)
Lagrange (11pt)| ~10^-10
BSpline (order 8)| ~10^-8
FastLagrange (11pt)| ~10^-6 (small shift ok)
```

## Choosing a Method

**Use Spectral or Lagrange if:**
- Maximum accuracy is required
- Function is very smooth (C^∞ or analytic)
- Small to medium grid sizes (< 10,000 points)

**Use BSpline if:**
- You want smooth interpolation with controlled properties
- Interpolating non-smooth data
- Need C² or C⁴ continuity

**Use FastLagrange if:**
- Performance is critical
- Shifts are small (< 1 grid spacing)
- Large-scale problems with many interpolations
- Memory limited

## Performance Tips

1. **Reuse objects**: Create interpolation objects once, use multiple times
   ```julia
   lagr = Lagrange(n, 5)
   for alpha in alphas
       interpolate!(f_interp, lagr, f, alpha)
   end
   ```

2. **Batch operations**: Process multiple arrays efficiently
   ```julia
   lagr = Lagrange(n, 5)
   for data in dataset
       interpolate!(output, lagr, data, shift)
   end
   ```

3. **Choose appropriate order**: Higher polynomial order = higher accuracy but slower
   - Order 3: Just for comparison, generally too low
   - Order 5: Good balance for most applications
   - Order 7+: Use only when high accuracy is critical

## Periodic Boundary Conditions

All methods assume **periodic boundary conditions**:
- Values are wrapped at grid boundaries
- `u[n+1] = u[1]`, `u[0] = u[n]`
- Large shifts automatically wrap around

```julia
n = 10
x = (0:n-1) ./ n
u = sin.(2π .* x)

lagr = Lagrange(n, 5)
u_out = zeros(n)

# These are equivalent due to periodicity:
interpolate!(u_out, lagr, u, 2.3)
interpolate!(u_out, lagr, u, 2.3 - 2.0)  # wraps to 0.3
```

## Error Analysis

For shift `α` and stencil of width `s` in the FastLagrange method:
- Truncation error ~ `O(h^p)` where `p = s-1` (polynomial order)
- Runge phenomena: weight magnitude grows as `~|α|^p` for `|α| > s/2`
- Best accuracy when `|α| < 0.5`

For spectral and FFT-based methods:
- Accuracy limited by function smoothness
- For C^k functions, error ~ `O(e^(-π k / n))` (exponential in smoothness)

## References

- Abramowitz, M., Stegun, I. A. (1964). Handbook of Mathematical Functions. Chapter 25.2
- FFTW Documentation: https://www.fftw.org/
- Julia FFTW.jl: https://github.com/JuliaComputing/AbstractFFTs.jl
