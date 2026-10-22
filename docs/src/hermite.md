# Hermite

The `Hermite` line interpolant accepts the same `bc` keyword as `Lagrange`:
`PeriodicBC()` is the default, while `ExtrapolationBC()` keeps stencils within
the data line and extrapolates at its ends.

```@autodocs
Modules = [SemiLagrangian]
Private = false
Order   = [:type, :function]
Pages   = ["hermite.jl", "hermite_interp.jl"]
```
