```@meta
CurrentModule = SemiLagrangian
```

# Spectral interpolation

[`Spectral`](@ref) shifts real-valued data on a uniform periodic grid using
Fourier transforms. It transforms the input to Fourier space, applies a phase
shift to each mode, and transforms the result back. It is well suited to smooth
periodic functions; discontinuities can produce Gibbs oscillations.

The current implementation is intended for an even number of grid points.
The periodic endpoint is excluded from the input. The shift is expressed in
grid cells: a positive shift samples the input to the right of each grid point.

## Example

This example shifts a smooth periodic function by a fraction of a grid cell and
compares the result with its exact value:

```@example spectral
using SemiLagrangian

n = 100
x = (0:n-1) ./ n
u = cos.(2pi .* x)
shift = 0.2
expected = cos.(2pi .* (x .+ shift / n))

interpolant = Spectral(n)
output = similar(u)
interpolate!(output, u, shift, interpolant)

println("Maximum error: ", maximum(abs.(output .- expected)))
```

Construct the interpolant once and reuse it for successive shifts. It owns
mutable FFT work arrays, so concurrent calls should use separate interpolants.

## Reference

```@docs
Spectral
interpolate!(::Any, ::Any, ::Real, ::Spectral)
```
