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
