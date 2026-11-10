```@meta
CurrentModule = SemiLagrangian
```

# FastInterpolations.jl

`FastInterpolations.jl` provides cubic-spline interpolation, including periodic
boundary conditions. It can be used directly with the uniform grids used in
the examples:

```julia
using FastInterpolations

n = 100
xi = 1.0:n
u = sin.(2π .* (0:n-1) ./ n)
alpha = 0.5

u_shifted = cubic_interp(
    xi,
    u,
    xi .+ alpha;
    bc=PeriodicBC(endpoint=:exclusive),
)
```

`endpoint=:exclusive` is appropriate for periodic grids that omit the duplicated
endpoint. The query coordinates are expressed in the same units as `xi`; here
`alpha` shifts the interpolated values by that many grid points. See
[`examples/landau-damping.jl`](https://github.com/JuliaVlasov/SemiLagrangian.jl/blob/master/examples/landau-damping.jl)
for a complete periodic Vlasov-Poisson example.

## CubicSpline wrapper

[`CubicSpline`](@ref) wraps `FastInterpolations.jl` with the package's
`interpolate!` interface for a uniform periodic grid. The input and output
arrays have `n` elements, with the duplicated periodic endpoint excluded.
The shift is measured in grid cells; a positive shift samples the input to
the right of each grid point.

```@example cubic_spline
using SemiLagrangian

n = 100
x = (0:n-1) ./ n
u = sin.(2pi .* x)
shift = 0.5
expected = sin.(2pi .* (x .+ shift / n))

interpolant = CubicSpline(n)
output = similar(u)
interpolate!(output, u, shift, interpolant)

println("Maximum error: ", maximum(abs.(output .- expected)))
```

### Reference

```@docs
CubicSpline
interpolate!(::Any, ::Any, ::Real, ::CubicSpline)
```
