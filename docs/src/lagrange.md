```@meta
CurrentModule = SemiLagrangian
```

# Lagrange

## Lagrange

`Lagrange(order)` constructs a general-purpose interpolant using `order + 1`
points per stencil.

`Lagrange` and `Hermite` accept a `bc` keyword. The default,
`PeriodicBC()`, wraps interpolation stencils around the ends of the data line.
Use `ExtrapolationBC()` to keep stencils within the line and extrapolate near
its ends:

```julia
using SemiLagrangian

periodic_interp = Lagrange(5; bc = PeriodicBC())
extrapolating_interp = Lagrange(5; bc = ExtrapolationBC())
```

```@docs
Lagrange
```

## LagrangeFFT

`LagrangeFFT(n, stencil)` applies a constant shift to a periodic line using
FFT-based convolution. It takes the grid length and the stencil width:

```julia
interpolant = LagrangeFFT(100, 7)
```

```@docs
LagrangeFFT
interpolate!(::Any, ::Any, ::Real, ::LagrangeFFT)
```

## LocalLagrange

`LocalLagrange(stencil)` applies a local fixed-width stencil to a periodic line.
The supported widths are 3, 5, 7, 9, and 11. It is intended for small shifts,
preferably within one grid cell.

```julia
interpolant = LocalLagrange(7)
```

```@docs
LocalLagrange
interpolate!(::Any, ::Any, ::Real, ::LocalLagrange)
```
