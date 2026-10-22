```@meta
CurrentModule = SemiLagrangian
```

# Precision

`SemiLagrangian` provides general-purpose line interpolants as well as
specialized methods for repeatedly shifting data on a one-dimensional,
uniform periodic grid. This page compares the specialized periodic methods
(`LagrangeFFT`, `Spectral`, `FastLagrange`, `CubicSpline`) with the
general-purpose interpolants (`Lagrange`, `Hermite`, `BSplineLU`,
`BSplineFFT`).

The specialized methods use the in-place interface:

```julia
interpolate!(output, interpolant, input, shift)
```

The general-purpose interpolants use the argument order
`interpolate!(output, input, shift, interpolant)`.

`shift` is measured in grid points, not physical-coordinate units. The
interpolated values correspond to sampling the input at `i + shift`; positive
shifts therefore move a sampled profile to the left. The input and output
must have the grid length used to construct the interpolant. Each method
treats the data as periodic, wrapping stencil indices or Fourier phases as
needed. Unlike the general `Lagrange` and `Hermite` interpolants, these
specialized types do not take a `bc` keyword: periodicity is built into them.

## Compare the methods

The following example compares each method with the exact shift of a smooth
periodic function. Each interpolant and output buffer is created once and
reused. The general-purpose interpolants are used with their default periodic
boundary condition (`PeriodicBC()`):

```julia
using SemiLagrangian

n = 128
x = (0:n-1) ./ n
u = sin.(2pi .* x) .+ 0.2 .* cos.(6pi .* x)
shift = 0.3
expected = sin.(2pi .* (x .+ shift / n)) .+
    0.2 .* cos.(6pi .* (x .+ shift / n))
output = similar(u)

methods = (
    "Lagrange FFT" => LagrangeFFT(n, 7),
    "Spectral" => Spectral(n),
    "Fast Lagrange" => FastLagrange(7),
    "Cubic spline" => CubicSpline(n),
)

for (name, interpolant) in methods
    interpolate!(output, interpolant, u, shift)
    println(rpad(name, 20), " max error = ", maximum(abs, output .- expected))
end

general = (
    "Lagrange (order 7)" => Lagrange(7),
    "Hermite (order 5)" => Hermite(5),
    "B-spline LU (order 7)" => BSplineLU(7, n),
    "B-spline FFT (order 7)" => BSplineFFT(7, n),
)

for (name, interpolant) in general
    interpolate!(output, u, shift, interpolant)
    println(rpad(name, 24), " max error = ", maximum(abs, output .- expected))
end
```

## Choosing a method

| Method | Approach | Consider it when |
| --- | --- | --- |
| `Spectral(n)` | Fourier phase shift | The data is smooth and high accuracy is important. Discontinuities can cause Gibbs oscillations. |
| `LagrangeFFT(n, stencil)` | Lagrange interpolation with FFT-based periodic convolution | You want a high-order periodic interpolant and can afford FFT workspaces. |
| `FastLagrange(stencil)` | Local, fixed-width polynomial stencil | Shifts are small and avoiding FFT setup is useful. |
| `CubicSpline(n)` | Cubic interpolation through the `FastInterpolations.jl` backend | A convenient periodic cubic spline is sufficient. |
| `Lagrange(order)` | Local Lagrange stencil, periodic or extrapolated boundaries | You need non-periodic boundaries or a different element type (`BigFloat`, `Rational`...). |
| `Hermite(order)` | Local Hermite interpolation using derivative estimates | You want a compact stencil with the Hermite construction (order 5, 9, ... or `flbis=true` with 3, 7, ...). |
| `BSplineLU(order, n)` | B-spline interpolation with an LU-factorised spline system | You need B-splines with periodic or extrapolated boundaries. |
| `BSplineFFT(order, n)` | Periodic B-spline interpolation solved with FFTs | You need periodic B-splines and want an FFT solve. |

The specialized interpolants allocate working arrays when constructed.
Reuse each object for successive time steps rather than constructing one per
call.

### Stencil and order choices

`LagrangeFFT` takes a grid length and a supported odd stencil width (for
example, `LagrangeFFT(n, 7)`). `FastLagrange` takes only the odd stencil
width; the implementation supports widths 3, 5, 7, 9, and 11. Its weights are
intended for small offsets, preferably within about one grid cell. For larger
shifts, split the time step or choose an FFT-based method.

## Reference

```@docs
LagrangeFFT
CubicSpline
Spectral
interpolate!(::Any, ::LagrangeFFT, ::Any, ::Real)
interpolate!(::Any, ::CubicSpline, ::Any, ::Real)
interpolate!(::Any, ::Spectral, ::Any, ::Float64)
```
