# Hermite

The `Hermite` line interpolant accepts the same `bc` keyword as `Lagrange`:
`PeriodicBC()` is the default, while `ExtrapolationBC()` keeps stencils within
the data line and extrapolates at its ends.

## Example

Shift a periodic function by a fraction of a grid cell and compare the result
with its exact value. The periodic grid excludes its repeated endpoint, and a
positive shift samples the input to the right of each grid point.

```@example hermite
using SemiLagrangian

n = 128
x = (0:n-1) ./ n
u = sin.(2pi .* x)
shift = 0.3
expected = sin.(2pi .* (x .+ shift / n))

interpolant = Hermite(5)
output = similar(u)
interpolate!(output, u, shift, interpolant)

println("Maximum error: ", maximum(abs.(output .- expected)))
```

## Reference

```@autodocs
Modules = [SemiLagrangian]
Private = false
Order   = [:type, :function]
Pages   = ["hermite.jl"]
```
